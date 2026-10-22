#' Helper function to get all accounting terms
#'
#' @param data tidyproteomics data object
#' @param format a character string specifying the return format: "table" (default) or "vector"
#'
#' @return a tibble or character vector
#' @export
#'
#' @examples
#' library(tidyproteomics)
#' get_accountings(hela_proteins)
#' get_accountings(hela_proteins, format = "vector")
#'
get_accountings <- function(
    data = NULL,
    format = c("table", "vector")
) {
  check_data(data)
  format <- rlang::arg_match(format)

  names_accounting <- names(data$accounting)
  names_accounting <- names_accounting[-which(grepl("sample|protein$|peptide$|modification$|group$", names_accounting))]

  if (format == "table") {
    data_types <- sapply(names_accounting, function(x) class(data$accounting[[x]])[1])
    tbl <- tibble::tibble(
      accounting = names_accounting,
      data_type = as.character(data_types)
    )
    return(tbl)
  }

  return(names_accounting)
}

