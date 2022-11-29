#' Helper function to take the head of a tibble and display as a data.frame
#'
#' @param x a tibble
#' @param n display up to the nth row
#'
#' @return a data frame
#' @export
#'
#' @examples
#' library(tidyproteomics)
#' x <- tibble::tibble(a = 1:10, b = 11:20)
#' hdf(x)
#' hdf(x, n = 3)
hdf <- function(x, n = 5) { as.data.frame(utils::head(x, n)) }
