#' Helper function to get all sample names
#'
#' @description
#' `get_sample_names()` is a helper function that returns the unique sample
#' names present in a tidyproteomics data object.
#'
#' @param data tidyproteomics data object
#'
#' @return a character vector of unique sample names
#' @export
#'
#' @examples
#' library(tidyproteomics)
#' get_sample_names(hela_proteins)
#'
get_sample_names <- function(
    data = NULL
) {
  check_data(data)
  data$experiments$sample %>% unique()
}
