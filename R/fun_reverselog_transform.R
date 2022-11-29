#' Reverse the plot axis for log transformation
#'
#' @param base a numeric
#'
#' @return a ggplot scale transformation
#'
reverselog_transformation <- function(base = exp(1)) {
  trans <- function(x) -log(x, base)
  inv <- function(x) base^(-x)
  scales::trans_new(paste0("reverselog-", format(base)), trans, inv,
                    scales::log_breaks(base = base),
                    domain = c(1e-100, Inf))
}
