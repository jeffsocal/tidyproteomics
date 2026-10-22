#' Helper function to get available terms
#'
#' @param data tidyproteomics data object
#' @param segment a character string
#' @param format a character string specifying the return format: "table" (default) or "vector"
#'
#' @return a tibble or character vector
#' @export
#'
#' @examples
#' library(tidyproteomics)
#' get_variables(hela_proteins, "experiments")
#' get_variables(hela_proteins, "experiments", format = "vector")
#'
get_variables <- function(
    data = NULL,
    segment = c('experiments', 'quantitative', 'annotations', 'accounting'),
    format = c("table", "vector")
) {
  check_data(data)
  format <- rlang::arg_match(format)

  segment <- rlang::arg_match(segment)
  cols_exp <- names(data$experiments)

  if(segment == 'experiments') {
    vars <- cols_exp
  } else if(segment == 'quantitative') {
    vars <- setdiff(names(data$quantitative), cols_exp)
  } else if(segment == 'accounting') {
    vars <- setdiff(names(data$accounting), c(data$identifier, cols_exp))
  } else if(segment == 'annotations') {
    vars <- unique(data$annotations$term)
  } else {
    vars <- c(data$identifier, unique(data$annotations$term))
  }

  if (format == "table") {
    return(tibble::tibble(
      segment = segment,
      variable = vars
    ))
  }

  return(vars)
}

#' Helper function to get available terms
#'
#' @param data tidyproteomics data object
#' @param variable a character string
#' @param .verbose a boolean
#'
#' @return a character
#'
get_segment <- function(
    data = NULL,
    variable = NULL,
    .verbose = TRUE
) {
  check_data(data)

  if(!is.character(variable)) {cli::cli_abort("`{variable}` must be a character")}

  for(segment in c('experiments', 'quantitative', 'annotations', 'accounting')){
    if(variable %in% get_variables(data, segment, format = 'vector')) {return(segment)}
  }
  if(.verbose == TRUE) {cli::cli_inform(c("!" = "Variable `{variable}` not found in data."))}
  return(NULL)
}

