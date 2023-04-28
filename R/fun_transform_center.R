#' helper function for normalizing a quantitative table
#'
#' @param table a tibble
#' @param group_by character vector
#' @param values character string
#' @param method character string
#'
#' @importFrom rlang :=
#'
#' @return a tibble
#'
center <- function(
    table,
    group_by = c("identifier"),
    values = "abundance",
    method = c('median','mean','geomean','sum')
){

  method <- rlang::arg_match(method)
  values <- rlang::arg_match(values, names(table))

  if(method == 'mean') {fun_center <- function(x) {mean(x, na.rm=T)}}
  if(method == 'median') {fun_center <- function(x) {stats::median(x, na.rm=T)}}
  if(method == 'geomean') {fun_center <- function(x) {exp(mean(log(x), na.rm=T))}}
  if(method == 'sum') {fun_center <- function(x) {sum(x, na.rm=T)}}

  value_method <- paste(values, method, sep="_")

  fun_app <- function(df, fun = mean, val = values, ...) {
    df %>% dplyr::summarise(dplyr::across(dplyr::starts_with(val), fun, !!!dplyr::enquos(...)), .groups = "drop")
  }

  table  %>%
    dplyr::group_by_at(group_by) %>%
    fun_app(fun = fun_center, val = values) %>%
    dplyr::rename({{value_method}} := {{values}})
}
