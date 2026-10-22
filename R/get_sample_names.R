#' Helper function to get all sample names
#'
#' @description
#' `get_sample_names()` is a helper function that returns the unique sample
#' names present in a tidyproteomics data object.
#'
#' @param data tidyproteomics data object
#' @param format a character string specifying the return format: "table" (default) or "vector"
#'
#' @return a tibble or character vector of unique sample names
#' @export
#'
#' @examples
#' library(tidyproteomics)
#' get_sample_names(hela_proteins)
#' get_sample_names(hela_proteins, format = "vector")
#'
get_sample_names <- function(
    data = NULL,
    format = c("table", "vector")
) {
  check_data(data)
  format <- rlang::arg_match(format)

  if (format == "table") {
    tbl <- data$experiments %>%
      dplyr::group_by(sample) %>%
      dplyr::summarise(
        replicates = dplyr::n_distinct(replicate),
        .groups = "drop"
      )
    return(tbl)
  }

  data$experiments$sample %>% unique()
}

