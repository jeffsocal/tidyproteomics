#' Imputes missing values based on a given function
#'
#' @param table a tibble
#' @param impute_function summary statistic function. Default is median, examples of
#' other functions include min, max, mean, sum. One could write a function for the
#' lower 5%, lower5th <- function(x) { quantile(x, 0.05)[[1]] }.
#' @param minimum_to_impute the minimum ratio to impute at
#'
#' @return a tibble
#'
impute_byfunction <- function(
    table = NULL,
    impute_function = stats::median,
    minimum_to_impute = 0.25
){

  # visible bindings

  identifier <- NULL
  n_groups <- NULL
  abundance <- NULL
  sample_rep <- NULL

  check_table(table)

  impute_function_name <- gsub("\\.", "_", as.character(methods::functionBody(impute_function))[2])

  tryCatch({

    tb_exps <- table %>%
      dplyr::select(sample, replicate) %>%
      unique() %>%
      dplyr::mutate(sample_rep = paste(sample, replicate, sep="_"))

    tb_keep <- table %>%
      dplyr::group_by(identifier) %>%
      dplyr::summarise(
        n_groups = length(identifier),
        .groups = 'drop') %>%
      dplyr::arrange(dplyr::desc(n_groups)) %>%
      dplyr::filter(n_groups >= tb_exps %>% nrow() * minimum_to_impute)

    table <- table %>% dplyr::filter(identifier %in% tb_keep$identifier)

    # pivot to wider table
    table_wide <- table %>%
      tidyr::pivot_wider(
        names_from = c('sample','replicate'),
        values_from = 'abundance',
        values_fill = NA
      ) %>%
      dplyr::mutate(row_id = dplyr::row_number() %>% as.character())

    # pivot back to collect/tag missing values
    table_long <- table_wide %>%
      tidyr::pivot_longer(
        cols = -c('identifier','row_id'),
        names_to = 'sample_rep',
        values_to = 'abundance'
      ) %>%
      dplyr::mutate(imputed = is.na(abundance))

    table_impute <- table_long %>%
      dplyr::group_by(sample_rep) %>%
      dplyr::mutate(abundance = ifelse(is.na(abundance),
                                       impute_function(abundance, na.rm=TRUE),
                                       abundance)) %>%
      dplyr::ungroup() %>%
      dplyr::inner_join(
        tb_exps, by='sample_rep'
      ) %>%
      dplyr::select("identifier", "sample", "replicate", "abundance", "imputed")

  }, error = function(err) {
    err = as.character(as.vector(err))
    cli::cli_abort(err)
  })

  return(table_impute)
}
