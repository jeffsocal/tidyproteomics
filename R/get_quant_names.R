#' Get the quantitative value names
#'
#' @description
#' `get_quant_names()` is a helper function that returns the names for all of the
#' normalized quantitative values, such as _raw_, _linear_, _loess_
#'
#' @param data a tidyproteomics data-object
#' @param format a character string specifying the return format: "table" (default) or "vector"
#'
#' @return a tibble or character vector
#' @export
#'
#' @examples
#' library(tidyproteomics)
#' get_quant_names(hela_proteins)
#' get_quant_names(hela_proteins, format = "vector")
#'
get_quant_names <- function(
    data = NULL,
    format = c("table", "vector")
) {
  check_data(data)
  format <- rlang::arg_match(format)

  abn_vals <- names(data$quantitative)[which(grepl('abundance', names(data$quantitative)))]
  abn_vals <- sub("abundance\\_", "", abn_vals)

  if (format == "table") {
    tbl <- tibble::tibble(
      quant_name = abn_vals,
      selected = (abn_vals == data$quantitative_source)
    )
    return(tbl)
  }

  return(abn_vals)
}

