#' Normalization function for a tidyproteomics data-object
#'
#' @param data tidyproteomics list data-object
#' @param data_centered a tibble of centered values used for normalization
#'
#' @return a tibble
#'
normalize_linear <- function(
    data = NULL,
    data_centered = NULL
){

  # visible bindings
  identifier <- NULL
  replicate <- NULL
  sample <- NULL
  abundance <- NULL
  abundance_centered <- NULL

  data_norm <- data %>%
    dplyr::left_join(data_centered, by='identifier', suffix = c("", "_centered")) %>%
    dplyr::select(identifier, replicate, sample, abundance, abundance_centered) %>%
    dplyr::group_by(sample, replicate) %>%
    tidyr::nest(data = c('identifier', 'abundance', 'abundance_centered')) %>%
    dplyr::mutate(m = purrr::map(data, lm, formula = abundance_centered ~ abundance + 0))

  for( i in 1:nrow(data_norm) ){
    tdf <- data_norm$data[[i]]
    lm <- data_norm$m[[i]]
    data_norm$data[[i]]$abundance_normalized <- stats::predict(lm, newdata=tdf)
  }

  nms <- c('identifier', 'sample', 'replicate', 'abundance_normalized')

  data_norm <- data_norm %>%
    tidyr::unnest(data) %>%
    dplyr::ungroup() %>%
    dplyr::select(nms)

  return(data_norm)
}
