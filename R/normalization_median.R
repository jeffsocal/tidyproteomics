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

  # compute the median shift
  data_centered <- data_centered %>%
    mutate(shift = median(abundance) - abundance)

  # apply the median shift and return the data
  data %>%
    left_join(data_centered,
              by = c('sample','replicate'),
              suffix = c("","_median")) %>%
    mutate(abundance_normalized = abundance + shift) %>%
    select(identifier, sample, replicate, abundance_normalized) %>%
    return()
}
