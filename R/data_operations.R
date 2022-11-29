#' Returns the data transformations
#'
#' @description
#' `operations()` returns the transformative operations performed on the data.
#'
#' @param data tidyproteomics data object
#' @param destination a character string
#'
#' @return a character
#' @export
#'
#' @examples
#' library(dplyr, warn.conflicts = FALSE)
#' library(tidyproteomics)
#' ecoli_proteins <- path_to_package_data("ecoli-hint_proteins") %>%
#'    import("ProteomeDiscoverer", "proteins") %>%
#'    reassign("sample", "hint_ko", "ko") %>%
#'    impute() %>%
#'    normalize(.method = c("linear","loess"))
#'
#' ecoli_proteins %>% operations()
#'
operations <- function(
    data=NULL,
    destination=c('print','save')
){

  # visible bindings
  n <- NULL

  destination <- rlang::arg_match(destination)
  check_data(data)
  v_trans <- data$operations %>% unlist() %>% as.vector()
  if(length(v_trans) == 0){
    cli::cli_alert_info("Data has not yet been transformed.")
  } else if(destination == "print") {
    cli::cli_inform(c("i" = "Data Transformations"))
    cli::cli_ol()
    ulid <- cli::cli_ul()
    cli::cli_li(v_trans)
    cli::cli_end(ulid)
    cli::cli_end()
  } else {
    tibble::tibble(data$operations %>% unlist() %>% as.vector()) %>%
      dplyr::mutate(n = dplyr::row_number()) %>%
      dplyr::relocate(n) %>%
      readr::write_csv("table_data_operations.csv")
  }
}
