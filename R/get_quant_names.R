#' Get the quantitative value names
#'
#' @description
#' `get_quant_names()` is a helper function that returns the names for all of the
#' normalized quantitative values, such as _raw_, _linear_, _loess_
#'
#' @param data a tidyproteomics data-object
#'
#' @return a character vector
#' @export
#'
#' @examples
#' library(tidyproteomics)
#' get_quant_names(ecoli_proteins)
#'
get_quant_names <- function(data) {
  check_data(data)
  abn_vals <- names(data$quantitative)[which(grepl('abundance', names(data$quantitative)))]
  abn_vals <- sub("abundance\\_", "", abn_vals)
  return(abn_vals)
}
