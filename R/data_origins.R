#' Assign sample origin to data experiments
#'
#' @description
#' `origins()` enables assigning or updating the `sample_origin` descriptive column
#' in the experimental table. This column is used by functions like `expression_me()`
#' (linear mixed-effects modeling) to designate biological replicate or subject origin
#' across runs. If `sample_origin` does not already exist, it is initialized to `sample_file`.
#'
#' @param data a tidyproteomics data-object
#' @param ... a three part expression (eg. sample_file == "F1", sample == "control")
#' @param .replace a character string
#'
#' @return a tidyproteomics data-object
#' @export
#'
#' @examples
#' library(dplyr, warn.conflicts = FALSE)
#' library(tidyproteomics)
#'
#' # check the experiment table
#' hela_proteins %>% summary("experiment")
#'
#' # assign sample origin (e.g. mapping files/runs to patient or donor IDs)
#' hela_proteins %>%
#'    origins(sample_file == "F1", .replace = "donor_1") %>%
#'    origins(sample_file == "F4", .replace = "donor_2") %>%
#'    summary("experiment")
#'
origins <- function(
    data = NULL,
    ...,
    .replace = NULL
){

  # visible bindings
  sample_id <- NULL

  check_data(data)

  # initialize sample_origin column if not present
  if(!'sample_origin' %in% colnames(data$experiments)){
    data$experiments$sample_origin <- data$experiments$sample_file
  }

  str_quo <- tidyproteomics_quo(...)
  if(is.null(str_quo)) { return(data) }

  if(is.null(.replace) || !is.character(.replace) || length(.replace) != 1) {
    cli::cli_abort(c("x" = "{.arg .replace} must be a single character string"))
  }

  variable <- str_quo[['variable']]
  operator <- str_quo[['operator']]
  value <- str_quo[['value']]
  inverse <- str_quo[['inverse']]
  inverse_str <- ''
  if(inverse == TRUE) { inverse_str <- '!' }

  which_segment <- get_segment(data, variable, .verbose = FALSE)
  if(is.null(which_segment) || which_segment != 'experiments'){
    cli::cli_div(theme = list(span.info = list(color = "#ff4500", "font-style" = 'italic')))
    cli::cli_abort(c("x" = "Variable {.info {variable}} not part of sample identifiers.\n Use one of {.emph {names(data$experiments)}}"))
  }

  which_rows <- down_select(data$experiments, str_quo)
  if(identical(which_rows, data$experiments)){
    # If down_select returned everything, verify if it truly matched all rows or matched none
    matches_all <- FALSE
    col_vals <- data$experiments[[variable]]
    if(operator == "==" && all(col_vals == value)) matches_all <- TRUE
    if(operator == "!=" && all(col_vals != value)) matches_all <- TRUE
    if(operator == "%like%" && all(grepl(value, col_vals, ignore.case = TRUE))) matches_all <- TRUE

    if(!matches_all){
      cli::cli_div(theme = list(span.emph = list(color = "#ff4500")))
      cli::cli_abort("did not find {.emph {value}} in {.emph {variable}}")
    }
  }

  w <- which(data$experiments$sample_id %in% which_rows$sample_id)
  if(length(w) == 0){
    cli::cli_div(theme = list(span.emph = list(color = "#ff4500")))
    cli::cli_abort("did not find matching rows for {.emph {variable}} {operator} {.emph {value}}")
  }

  cli::cli_alert_info('Setting origin for {length(w)} run{?s}')

  data$experiments$sample_origin[w] <- as.character(.replace)

  data$operations <- append(data$operations, glue::glue("Data assigned sample_origin to '{.replace}' where '{variable} {operator} {value}'"))
  return(data)
}
