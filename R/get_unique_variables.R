#' Helper function to get all sample names
#'
#' @param data tidyproteomics data object
#' @param variable a string character
#'
#' @return a vector
#'
get_unique_variables <- function(
    data = NULL,
    variable = NULL
) {

  # visible bindings
  term <- NULL
  annotation <- NULL

  check_data(data)

  which_segment <- get_segment(data, variable)
  if(!is.null(which_segment)) {

    if(which_segment == 'annotations') {
      data[[which_segment]] <- data[[which_segment]] %>%
        dplyr::filter(term == variable) %>%
        tidyr::separate_rows(annotation, sep=";")
      variable <- 'annotation'
    }
      unique(unlist(data[[which_segment]][,variable]))
  }

}
