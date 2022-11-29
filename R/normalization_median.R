#' Normalization function for a tidyproteomics data-object
#'
#' @param data tidyproteomics list data-object
#' @param data_centered a tibble of centered values used for normalization
#'
#' @return a tibble
#'
normalize_median <- function(
    data = NULL,
    data_centered = NULL
){

  # visible bindings
  abundance <- NULL

  v_med <- data_centered$abundance %>% stats::median()

  data_norm <- list()

  for( i in 1:nrow(data_centered) ){
    tdf <- data %>% dplyr::filter(sample == data_centered$sample[i] & replicate == data_centered$replicate[i])
    t_med <- data_centered$abundance[i]
    data_norm[[i]] <- tdf %>%
      dplyr::mutate(abundance_normalized = abundance + (v_med - t_med))
  }

  nms <- c('identifier', 'sample', 'replicate', 'abundance_normalized')

  data_norm <- data_norm %>%
    dplyr::bind_rows() %>%
    dplyr::select(nms)

  return(data_norm)
}
