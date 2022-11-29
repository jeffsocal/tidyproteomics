#' Calculates the sum of a numeric vector with NAs removed
#'
#' @param x a numeric vector
#' @return a numeric
#' @export
#'
#' @examples
#' library(tidyproteomics)
#' fsum(c(1,2,5,6,8,NA,NA))
#'
fsum <- function(x) {sum(x, na.rm=T)}

#' Calculates the mean of a numeric vector with NAs removed
#'
#' @param x a numeric vector
#' @return a numeric
#' @export
#'
#' @examples
#' library(tidyproteomics)
#' fmean(c(1,2,5,6,8,NA,NA))
#'
fmean <- function(x) {mean(x, na.rm=T)}

#' Calculates the median of a numeric vector with NAs removed
#'
#' @param x a numeric vector
#' @return a numeric
#' @export
#'
#' @examples
#' library(tidyproteomics)
#' fmedian(c(1,2,5,6,8,NA,NA))
#'
fmedian <- function(x) {stats::median(x, na.rm=T)}

#' Calculates the geometric mean of a numeric vector with NAs removed
#'
#' @param x a numeric vector
#' @return a numeric
#' @export
#'
#' @examples
#' library(tidyproteomics)
#' fgeomean(c(1,2,5,6,8,NA,NA))
#'
fgeomean <- function(x) {exp(mean(log(x[which(!is.na(x))])))}

#' Calculates the minimum of a numeric vector with NAs removed
#'
#' @param x a numeric vector
#' @return a numeric
#' @export
#'
#' @examples
#' library(tidyproteomics)
#' fmin(c(1,2,5,6,8,NA,NA))
#'
fmin <- function(x) {base::min(x, na.rm=T)}
