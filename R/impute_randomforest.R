#' Imputes missing values based on the missForest function
#'
#' @param table a tibble
#' @param minimum_to_impute the minimum ratio to impute at
#' @param cores the number of threads used to speed the calculation
#'
#' @return a tibble
#'
impute_randomforest <- function(
    table = NULL,
    minimum_to_impute = 0.25,
    cores = 2
){

  # visible bindings
  identifier <- NULL
  n_groups <- NULL
  abundance <- NULL
  sample_rep <- NULL
  imputed <- NULL
  dt_long <- NULL
  row_id <- NULL
  oobe <- NULL
  n_imputed <- NULL
  observation <- NULL

  # check_table(table)

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


    doParallel::registerDoParallel(cores=cores)
    # impute the missing values (ie. the NAs)
    object_impute <- table_wide %>%
      dplyr::select(-c('identifier')) %>%
      as.data.frame() %>%
      tibble::column_to_rownames('row_id') %>%
      t() %>%
      as.matrix() %>%
      missForest::missForest(
        parallelize = 'variables',
        variablewise = TRUE)

    table_oobe <- tibble::tibble(
      oobe = object_impute$OOBerror) %>%
      dplyr::mutate(row_id = dplyr::row_number() %>% as.character()) %>%
      dplyr::inner_join(
        table_long %>% dplyr::select(identifier, row_id, imputed),
        by = c('row_id')
      ) %>%
      dplyr::group_by(identifier) %>%
      dplyr::summarise(
        oobe = mean(oobe),
        n_imputed = length(which(imputed == T)),
        r_imputed = n_imputed / length(imputed),
        .groups = 'drop'
      )

    table_impute <- object_impute$ximp %>%
      t() %>%
      as.data.frame() %>%
      tibble::rownames_to_column('row_id') %>%
      tibble::as_tibble() %>%
      tidyr::pivot_longer(
        cols = !dplyr::matches('row_id'),
        names_to = 'sample_rep',
        values_to = 'abundance'
      ) %>%
      dplyr::inner_join(
        table_long %>% dplyr::select(!abundance),
        by = c('row_id','sample_rep')
      ) %>%
      dplyr::inner_join(
        tb_exps, by='sample_rep'
      ) %>%
      dplyr::select("identifier", "sample", "replicate", "abundance", "imputed")

  }, error = function(err) {
    cli::cli_abort(c("x" = as.character(as.vector(err))))
  })

  msg_echo <- table_oobe %>% dplyr::filter(n_imputed > 0) %>%
    dplyr::summarize(n = dplyr::n(), oobe = stats::median(oobe))

  l_out <- list(
    table = table_impute,
    operation = glue::glue('... median oobe {signif(msg_echo$oobe[[1]], 3)}')
  )

  cli::cli_alert_info(c('i' = l_out$operation))

  return(l_out)
}
