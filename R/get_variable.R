#' Helper function to get available terms
#'
#' @param data tidyproteomics data object
#' @param segment a character string
#'
#' @return a vector
#'
get_variables <- function(
    data = NULL,
    segment = c('experiments', 'quantitative', 'annotations', 'accounting')
) {
  check_data(data)

  segment <- rlang::arg_match(segment)
  cols_exp <- names(data$experiments)

  if(segment == 'experiments') {return(cols_exp)}
  if(segment == 'quantitative') {return(setdiff(names(data$quantitative), cols_exp))}
  if(segment == 'accounting') {return(setdiff(names(data$accounting), c(data$identifier, cols_exp)))}
  if(segment == 'annotations') {return(unique(data$annotations$term))}

  return(c(data$identifier, unique(data$annotations$term)))
}

#' Helper function to get available terms
#'
#' @param data tidyproteomics data object
#' @param variable a character string
#' @param verbose a boolean
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
    if(variable %in% get_variables(data, segment)) {return(segment)}
  }
  if(.verbose == TRUE) {cli::cli_inform(c("!" = "Variable `{variable}` not found in data."))}
  return(NULL)
}
