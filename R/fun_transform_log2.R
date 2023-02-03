#' helper function for normalizing a quantitative table
#'
#' @param table a tibble
#' @param values a character string
#'
#' @importFrom rlang :=
#'
#' @return a tibble
#'
transform_log2 <- function(
    table,
    values = "abundance"
    ){

  # visible bindings
  transform_values <- NULL

  values <- rlang::arg_match(values, names(table))

  colnames(table)[which(colnames(table) == values)] <- 'transform_values'

  values <- paste0(values, "_log2")

  table %>%
    dplyr::mutate(transform_values = ifelse(transform_values == 0, NA, transform_values)) %>%
    dplyr::filter(!is.na(transform_values)) %>%
    dplyr::mutate({{values}} := log2(transform_values)) %>%
    dplyr::select(!transform_values)
}
