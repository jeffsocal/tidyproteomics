#' Read data by format type
#'
#' @description
#' `read_data()` is a helper function that assumes the format type of the data
#' table by checking the ending of path string
#'
#' @param path a path character string
#' @param platform a character string
#' @param analyte a character string
#'
#' @return tibble
#'
read_data <- function(
    path = NULL,
    platform = NULL,
    analyte = c("peptides","proteins")
) {

  # visible bindings
  col_name <- NULL

  analyte <- rlang::arg_match(analyte)
  format <- c("txt", "csv", "tsv", "rds", "xlsx", "xls", "mzTab")
  ext <- stringr::str_extract(path, "\\..{3,6}$") %>% stringr::str_remove("^\\.")
  ext <- rlang::arg_match(ext, format)

  cli::cli_alert_info("... reading data as `{ext}`")

  if( ext == 'rds' ) {
    tbl <- path %>% readr::read_rds()
  } else if( ext == 'mzTab' && platform == 'mzTab' ) {
    obj <- path %>% read_mzTab(analyte)
    analyte <- obj$analyte
    platform <- obj$platform
    tbl <- obj$data
  } else if( ext %in% format[1:3] ) {
    tbl <- path %>% readr::read_tsv()
    cli::cli_alert_info('... data dimentions `{dim(tbl)}`')

    tbl_problems <- vroom::problems(tbl)
    if(nrow(tbl_problems) > 1){

      tbl_problems <- tbl_problems %>%
        dplyr::inner_join(
          tibble::tibble(col_name = colnames(tbl)) %>%
            dplyr::mutate(col = dplyr::row_number()),
          by='col') %>%
        dplyr::select(!dplyr::matches('file|col$')) %>%
        dplyr::relocate(col_name)

      cli::cli_div(theme = list(span.emph = list(color = "#ff4500")))
      cli::cli_alert_info("... issues with {.emph {signif(length(unique(tbl_problems$row))/nrow(tbl) * 100,3)}%} of the data imported")
      cli::cli_alert_info("... effecting columns {.emph {unique(tbl_problems$col_name)}}")
      cli::cli_alert_info("... check to verify that the columns effected are not imported")

      cli::cli_alert_info("{.emph ===== EXAMPLE =====}")
      print(tbl_problems %>% utils::head())
      cli::cli_alert_info("{.emph ===================}")
    }

  } else if( ext %in% format[5:6] ) {
    tbl <- path %>% readxl::read_excel(col_types = 'text')
  } else {
    cli::cli_div(theme = list(span.emph = list(color = "#ff4500")))
    cli::cli_abort("Can not import {.emph {path}}. File format not recognized, must be one of {.emph {format}}.")
  }

  return(list(platform = platform, analyte = analyte, data = tbl))
}

#' Load project specific data
#'
#' @description
#' `load_local()` is a simple function that loads the current project
#' tidyproteomics data object
#'
#' @param analyte a character string
#'
#' @return an tidyproteomics data object
#' @export
#'
#' @examples
#' library(tidyproteomics)
#' # hela_proteins <- load_omics(analyte = "proteins")
#'
load_local <- function(
    analyte = c("peptides","proteins")
) {
  analyte <- rlang::arg_match(analyte)
  path <- glue::glue("./data/rds/{analyte}.rds")
  if(!file.exists(path)) {
    cli::cli_abort(c("x" = "Data not found for {analyte}"))
  }
  return(path %>% readRDS())
}
