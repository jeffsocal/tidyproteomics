#' A function for evaluating expression differences between two sample sets via the limma algorithm
#'
#' @param data tidyproteomics data object
#' @param impute_max a numeric representing the largest allowable imputation percentage
#'
#' @return a tibble
#'
analysis_counts <- function(
    data = NULL,
    impute_max = 0.5
){

  # visible bindings
  match_between_runs <- NULL
  imputed <- NULL
  quantifiable <- NULL

  check_data(data)
  if(!is.numeric(impute_max)) {cli::cli_abort("{.emph impute_max} is not a numeric ..")}
  tbl_experiments <- data$experiments
  limit <- nrow(tbl_experiments)
  get_cols_experiments <- intersect(colnames(tbl_experiments), c('sample_id', 'sample', 'replicate', 'sample_group'))

  fat <- list(
    data %>%
      subset(match_between_runs == FALSE, .verbose = FALSE) %>%
      subset(imputed <= !!impute_max, .verbose = FALSE) %>%
      summary('sample_id', destination = 'return', limit = limit) %>% dplyr::mutate(is_mbr = FALSE),
    data %>% summary('sample_id', destination = 'return', limit = limit) %>% dplyr::mutate(is_mbr = TRUE)
  )  %>%
    dplyr::bind_rows() %>%
    dplyr::full_join(
      tbl_experiments %>% dplyr::select(get_cols_experiments),
      by = 'sample_id'
    ) %>%
    dplyr::mutate(quantifiable = quantifiable/100) %>%
    dplyr::relocate(get_cols_experiments) %>%
    dplyr::select(!dplyr::matches('CVs|files'))

  return(fat)
}
