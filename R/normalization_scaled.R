#' Normalization function for a tidyproteomics data-object
#'
#' @param data tidyproteomics list data-object
#' @param data_centered a tibble of centered values used for normalization
#'
#' @return a tibble
#'
normalize_scaled <- function(
    data = NULL,
    data_centered = NULL
){

  # visible bindings
  abundance <- NULL
  shift <- NULL
  identifier <- NULL
  abundance_normalized <- NULL

  # compute the scaled shift
  data_centered <- data_centered %>%
    dplyr::mutate(shift = stats::median(abundance, na.rm = TRUE) - abundance)

  # apply the shift and return the data
  data %>%
    dplyr::left_join(data_centered,
                     by = c('sample', 'replicate'),
                     suffix = c("", "_scaled")) %>%
    dplyr::mutate(abundance_normalized = abundance + shift) %>%
    dplyr::select(identifier, sample, replicate, abundance_normalized)
}
