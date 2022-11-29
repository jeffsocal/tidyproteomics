#' Helper function to get all sample names
#'
#' @param data tidyproteomics data object
#'
#' @return a vector
#'
get_sample_names <- function(
    data = NULL
) {
  check_data(data)
  data$experiments$sample %>% unique()
}
