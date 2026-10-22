#' Helper function to get unique values for a variable
#'
#' @param data tidyproteomics data object
#' @param variable a string character
#' @param format a character string specifying the return format: "table" (default) or "vector"
#'
#' @return a tibble or vector
#' @export
#'
#' @examples
#' library(tidyproteomics)
#' get_unique_variables(hela_proteins, "sample")
#' get_unique_variables(hela_proteins, "sample", format = "vector")
#'
get_unique_variables <- function(
    data = NULL,
    variable = NULL,
    format = c("table", "vector")
) {

  # visible bindings
  term <- NULL
  annotation <- NULL

  check_data(data)
  format <- rlang::arg_match(format)

  which_segment <- get_segment(data, variable, .verbose = FALSE)
  if(!is.null(which_segment)) {

    if(which_segment == 'annotations') {
      data[[which_segment]] <- data[[which_segment]] %>%
        dplyr::filter(term == variable) %>%
        tidyr::separate_rows(annotation, sep=";")
      variable <- 'annotation'
    }
    vals <- unique(unlist(data[[which_segment]][,variable]))

    if (format == "table") {
      tbl <- tibble::tibble(value = vals)
      return(tbl)
    }

    return(vals)
  }

  if (format == "table") {
    return(tibble::tibble(value = character()))
  }
  return(NULL)
}

