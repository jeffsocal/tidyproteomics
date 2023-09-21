#' Helper function to get all annotations for a given term
#'
#' @param data tidyproteomics data object
#' @param term a character string
#'
#' @return a vector
#'
get_annotations <- function(
    data = NULL,
    term = NULL
) {

  # visible bindings
  term_group <- term
  annotation <- NULL

  check_data(data)
  data$annotations %>%
    dplyr::filter(term == term_group) %>%
    tidyr::separate_rows(annotation, sep = ";") %>%
    dplyr::mutate(annotation = trimws(annotation)) %>%
    dplyr::select(!term)
}

#' Helper function to get available terms
#'
#' @param data tidyproteomics data object
#'
#' @return a vector
#'
get_annotation_terms <- function(data) {
  check_data(data)
  return(c(data$identifier, unique(data$annotations$term)))
}
