#' Analysis tables and plots of expression values
#'
#' @description
#' `analyze_expressions()` is a GGplot2 implementation for plotting the expression differences
#' as foldchange ~ statistical significance. See also `plot_proportion()`.  This function can
#' take either a tidyproteomics data object or a table with the required headers.
#'
#' @param log2fc_column a character defining the column name of the log2 foldchange values.
#' @param log2fc_min a numeric defining the minimum log2 foldchange to highlight.
#' @param significance_column a character defining the column name of the statistical significance values.
#' @param significance_max a numeric defining the maximum statistical significance to highlight.
#' @param labels_column a character defining the column name of the column for labeling.
#' @param show_pannels a boolean for showing colored up/down expression panels.
#' @param show_lines a boolean for showing threshold lines.
#' @param show_fc_scale a boolean for showing the secondary foldchange scale.
#' @param point_size a character reference to a numerical value in the expression table
#' @param color_positive a character defining the color for positive (up) expression.
#' @param color_negative a character defining the color for negative (down) expression.
#' @param height a numeric
#' @param width a numeric

#'
#' @return a tidyproteomics data object
#' @export
#'
#' @examples
#' library(dplyr, warn.conflicts = FALSE)
#' library(tidyproteomics)
#' hela_proteins %>%
#'    expression(knockdown/control) %>%
#'    analyze_expressions(log2fc_min = 0.5, significance_column = "p_value")
#'
analyze_expressions <- function(
    data = NULL,
    log2fc_min = 1,
    log2fc_column = 'log2_foldchange',
    significance_max = 0.05,
    significance_column = 'adj_p_value',
    labels_column = NULL,
    show_pannels = TRUE,
    show_lines = TRUE,
    show_fc_scale = TRUE,
    point_size = NULL,
    color_positive = 'dodgerblue',
    color_negative = 'firebrick1',
    height = 5,
    width = 8
){

  # visible bindings
  metric <- NULL

  check_data(data)

  # pull annotations if they exist
  if(!is.null(data$annotations)){
    tbl_anno <- data$annotations %>% tidyr::pivot_wider(names_from = 'term', values_from = 'annotation')
  }

  for( set_expression in names(data$analysis)){

    if(is.null(data$analysis[[set_expression]]$expression)){ next }

    table <- data$analysis[[set_expression]]$expression %>%
      dplyr::left_join(
        extract(data) %>%
          tidyr::pivot_wider(names_from = c('sample','replicate'), values_from = 'abundance', names_prefix = 'abundance_') %>%
          dplyr::rename(normalization = origin) %>%
          munge_identifier("separate", data$identifier),
        by = data$identifier)


    if(!is.null(data$annotations))
      table <- table %>% dplyr::left_join(tbl_anno, by = 'protein')

    table |> plot_volcano(
      log2fc_min = log2fc_min,
      log2fc_column = log2fc_column,
      significance_max = significance_max,
      significance_column = significance_column,
      labels_column = labels_column,
      show_pannels = show_pannels,
      show_lines = show_lines,
      show_fc_scale = show_fc_scale,
      point_size = point_size,
      color_positive = color_positive,
      color_negative = color_negative,
      destination = 'png',
      height = height,
      width = width
    )

    # plot the expressions
    plot_name <- glue::glue("{data$analyte}_volcano_{stringr::str_replace_all(set_expression, '/', '-')}.png")
    file.rename('plot_volcano.png', plot_name)

    # save the expression table
    data_name <- glue::glue("table_{data$analyte}_expression_{stringr::str_replace_all(set_expression, '/', '-')}.csv")
    table |> readr::write_csv(data_name)

  }

  return(data)
}
