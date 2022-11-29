#' Check the integrity of a tidyproteomics quantitative tibble
#'
#' @description
#' `check_table()` is a helper function that checks the structure and contents of
#' a tidyproteomics quantitative tibble
#'
#' @param table a tibble
#'
#' @return silent on success, an abort message on fail
#'
check_table <- function(
    table = NULL){

  # fail if data is NULL
  if(is.null(table)) {
    cli::cli_div(theme = list(span.emph = list(color = "#ff4500"), span.info = list(color = "red")))
    cli::cli_abort(c("x" = "Input {.emph table-object} not recognized", "{.info ... analysis terminated}"))
  }
  names_need <- c("identifier","sample","replicate","abundance")
  names_have <- names(table)
  names_inter <- intersect(names_need, names_have)
  if(length(names_inter) != 4) {
    # fail if missing table objects
    names_diff <- paste(dplyr::setdiff(names_need, names_have), collapse = ", ")
    cli::cli_div(theme = list(span.emph = list(color = "#ff4500")))
    cli::cli_alert_danger("Input {.emph table-object} not recognized")
    cli::cli_div(theme = list(span.emph = list(color = "red"), span.info = list(color = "red")))
    cli::cli_abort(c("x" = "Missing ({.emph {names_diff}}) in table-object", "{.info ... analysis terminated}"))
  }

  if(length(which(grepl("^[0-9]", table$sample, perl = TRUE))) > 0){
    cli::cli_div(theme = list(span.emph = list(color = "#ff4500"), span.info = list(color = "red")))
    cli::cli_abort(c("x" = "Input {.emph sample} names can not start with a  {.emph numeric}", "{.info ... analysis terminated}"))
  }

  if(length(which(grepl("\\-", table$sample, perl = TRUE))) > 0){
    cli::cli_div(theme = list(span.emph = list(color = "#ff4500"), span.info = list(color = "red")))
    cli::cli_abort(c("x" = "Input {.emph sample} names must not contain a  {.emph hyphen (-)}", "{.info ... analysis terminated}"))
  }
}
