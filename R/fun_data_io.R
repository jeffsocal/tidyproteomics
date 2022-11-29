#' Read data by format type
#'
#' @description
#' `read_data()` is a helper function that assumes the format type of the data
#' table by checking the ending of path string
#'
#' @param path a path character string
#' @param ... readr passthrough options
#'
#' @return tibble
#'
read_data <- function(
    path = NULL,
    ...
) {

  ext <- stringr::str_extract(path, "\\..{3,4}$")
  ext <- rlang::arg_match(ext, c(".rds", ".csv", ".tsv", ".xlsx", ".xls", ".txt"))

  if( grepl("\\.rds$", path) ) {
    path %>% readr::read_rds()
  } else if( grepl("\\.csv$", path) ) {
    path %>% readr::read_csv(...)
  } else if( grepl("\\.(tsv|txt)$", path) ) {
    path %>% readr::read_tsv(...)
  } else if( grepl("\\.xl[xs]+$", path) ) {
    path %>% readxl::read_excel()
  }
}

#' Load project specific data
#'
#' @description
#' `load_local()` is a simple function that loads the current project
#' tidyproteomics data object
#'
#' @param source a character string
#'
#' @return an tidyproteomics data object
#' @export
#'
#' @examples
#' library(tidyproteomics)
#' # ecoli_proteins <- load_omics(source = "proteins")
#'
load_local <- function(
    source = c("peptides","proteins")
) {
  source <- rlang::arg_match(source)
  path <- glue::glue("./data/rds/{source}.rds")
  if(!file.exists(path)) {
    cli::cli_abort(c("x" = "Data not found for {source}"))
  }
  return(path %>% readRDS())
}
