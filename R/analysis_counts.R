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
  limit <- nrow(data$experiments)

  fat <- list(
    data %>%
      subset(match_between_runs == FALSE, .verbose = FALSE) %>%
      subset(imputed <= !!impute_max, .verbose = FALSE) %>%
      summary('sample_id', destination = 'return', limit = limit) %>% dplyr::mutate(is_mbr = FALSE),
    data %>% summary('sample_id', destination = 'return', limit = limit) %>% dplyr::mutate(is_mbr = TRUE)
  )  %>%
    dplyr::bind_rows() %>%
    dplyr::full_join(
      data$experiment %>%
        dplyr::select(c('sample_id', 'sample', 'replicate')),
      by = 'sample_id'
    ) %>%
    dplyr::mutate(quantifiable = quantifiable/100) %>%
    dplyr::relocate('sample', 'replicate', 'is_mbr') %>%
    dplyr::select(!dplyr::matches('CVs|files'))

  return(fat)
}
