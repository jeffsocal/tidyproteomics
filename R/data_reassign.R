#' reassign the sample info
#'
#' @description
#' `reassign()` enables editing of the sample descriptive in the experimental table.
#' This function will only replace the sample string and update the replicate number.
#'
#' @param data a tidyproteomics data-object
#' @param ... a three part expression (eg. x == a)
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
#' # make the modification
#' hela_proteins %>%
#'    reassign(sample == "control", .replace = "ct") %>%
#'    reassign(sample == "knockdown", .replace = "kd") %>%
#'    summary("sample")
#'
#' # reassign specific file_ids
#' hela_proteins %>%
#'    reassign(sample_file == "f1", .replace = "new") %>%
#'    reassign(sample_file == "f2", .replace = "new") %>%
#'    summary("sample")
#'
reassign <- function(
    data = NULL,
    ...,
    .replace = NULL
){

  # visible bindings
  sample_new <- NULL
  replicate_new <- NULL

  # visible bindings
  sample_id <- NULL

  check_data(data)
  str_quo <- tidyproteomics_quo(...)
  if(is.null(str_quo)) { return(data) }

  variable <-  str_quo[['variable']]
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

  if(grepl("^[0-9]", .replace, perl = TRUE)){
    cli::cli_div(theme = list(span.emph = list(color = "#ff4500")))
    cli::cli_abort(c("x" = "Input {.emph sample} names must not start with a {.emph numeric}"))
  }

  if(grepl("/", .replace, perl = TRUE)){
    cli::cli_div(theme = list(span.emph = list(color = "#ff4500")))
    cli::cli_abort(c("x" = "Input {.emph sample} names must not contain a {.emph /}"))
  }

  which_rows <- down_select(data$experiments, str_quo)
  if(identical(which_rows, data$experiments)){
    cli::cli_div(theme = list(span.emph = list(color = "#ff4500")))
    cli::cli_abort("did not find {.emph {value}} in {.emph {variable}}")
  }

  w <- which(data$experiments$sample_id %in% which_rows$sample_id)

  cli::cli_alert_info('Reassigning {w}')

  data$experiments$sample[w] <- .replace
  data$experiments <- data$experiments %>%
    dplyr::group_by(sample) %>%
    dplyr::mutate(replicate = as.character(dplyr::row_number())) %>%
    dplyr::ungroup()

  data$quantitative <- data$quantitative %>%
    dplyr::select(!c("sample","replicate")) %>%
    dplyr::full_join(
      data$experiments %>% dplyr::select("sample_id","sample","replicate"),
      by = "sample_id"
    ) %>%
    dplyr::relocate(c('sample_id', 'sample', 'replicate', data$identifier))

  data$operations <- append(data$operations, glue::glue("Data reassigned sample to '{.replace}' where '{variable} {operator} {value}'"))
  return(data)
}
