#' Bubble plot of enrichment values
#'
#' @description
#' `plot_enrichment()` is a GGplot2 implementation for plotting the enrichment values. This function can
#' take either a tidyproteomics data object or a table with the required headers.
#'
#' @param data a tidyproteomics data object
#' @param ... two sample comparison
#' @param .term a character string indicating the term enrichment analysis should be calculated for
#' @param enrichment_column a character defining the column name of enrichment values.
#' @param enrichment_min a numeric defining the minimum log2 enrichment to highlight.
#' @param significance_column a character defining the column name of the statistical significance values.
#' @param significance_max a numeric defining the maximum statistical significance to highlight.
#' @param term_column a character defining the column name for labeling.
#' @param size_column a character defining the column name of term size.
#' @param destination a character string
#' @param height a numeric
#' @param width a numeric
#'
#' @return a ggplot2 object
#' @export
#'
#' @examples
#' library(dplyr, warn.conflicts = FALSE)
#' library(ggplot2, warn.conflicts = FALSE)
#' library(tidyproteomics)
#' hela_proteins %>%
#'    expression(knockdown/control, .method = stats::t.test) %>%
#'    enrichment(knockdown/control, .terms = 'biological_process', .method = "wilcoxon") %>%
#'    plot_enrichment(knockdown/control, .term = "biological_process") +
#'    labs(title = "Hela: Term Enrichment", subtitle = "Knockdown ~ Control")
#'
#'
plot_enrichment <- function(
    data = NULL,
    ...,
    .term = NULL,
    enrichment_min = 1,
    enrichment_column = 'enrichment',
    significance_max = 0.01,
    significance_column = 'p_value',
    term_column = 'annotation',
    size_column = 'size',
    destination = 'plot',
    height = 5,
    width = 8
){

  # visible bindings
  .data <- NULL
  term_col <- NULL
  keep <- NULL
  plot_title <- ""
  plot_subtitle <- ""

  file_name = 'enrichment'
  if('tidyproteomics' %in% class(data)) {
    table <- data %>% export_analysis(..., .analysis = 'enrichment', .term = .term)
    str_quo <- tidyproteomics_quo_name(...)
    file_name = glue::glue("{data$analyte}_{data$quantitative_source}_enrichment_{str_quo}")
    plot_title <- tidyproteomics_quo_name(..., sep = " / ") %>% stringr::str_to_title()
    plot_subtitle <- data$analysis[[tidyproteomics_quo_name(..., sep = "/")]]$enrichment[[.term]]$method
  } else {
    table <- data
  }

  table_cols <- colnames(table)
  enrichment_column <- rlang::arg_match(enrichment_column, table_cols)
  significance_column <- rlang::arg_match(significance_column, table_cols)
  term_column <- rlang::arg_match(term_column, table_cols)
  size_column <- rlang::arg_match(size_column, table_cols)

  table <- table %>%
    dplyr::filter(!is.na(.data[[enrichment_column]]),
                  !is.na(.data[[significance_column]])) %>%
    dplyr::mutate(keep = abs(.data[[enrichment_column]]) >= enrichment_min &
                    .data[[significance_column]] <= significance_max,
                  term_col = .data[[term_column]]) %>%
    dplyr::arrange(dplyr::desc(.data[[enrichment_column]])) %>%
    dplyr::mutate(term_col = forcats::fct_reorder(term_col, .data[[enrichment_column]], .desc = FALSE))


  # construct the enrichment plot
  plot <- table %>%
    ggplot2::ggplot(ggplot2::aes(x = .data[[enrichment_column]],
                                 y = term_col,
                                 color = .data[[significance_column]],
                                 size = .data[[size_column]])) +
    ggplot2::geom_vline(xintercept = 0, alpha=.33) +
    ggplot2::geom_point(alpha = 0.2) +
    ggplot2::geom_point(data = table %>% dplyr::filter(keep == TRUE),
                        alpha = 1) +
    ggplot2::scale_color_gradient(low = 'red', high = "blue") +
    ggplot2::theme_bw() +
    ggplot2::labs(x = "Enrichment",
                  y = "",
                  title = plot_title,
                  subtitle = plot_subtitle)

  return(plot_save(plot,
                   data,
                   file_name = file_name,
                   destination = destination,
                   height = height,
                   width = width))
}
