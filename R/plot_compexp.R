#' Comparative analysis between two expression tests
#'
#' @description
#' `plot_compexp()` is a GGplot2 implementation for plotting the comparison in
#' expression differences between two methods or two sets of groups. For example,
#' one could run an expression difference for two different conditions (A and B)
#' prodived the experiment contained 3 samples condition A, condition B and WT,
#' then compare those results. The proteins showing up in the intersection (purple)
#' indicate common targets for condition A and B.
#' ```
#' expdiff_a <- protein_data %>%
#'    expression(experiment = "condition_a", control = "wt")
#'
#' expdiff_b <- protein_data %>%
#'    expression(experiment = "condition_b", control = "wt")
#'
#' plot_compexp(expdiff_a, expdiff_b)
#'
#' ```
#'
#' @param table_a a tibble
#' @param table_b a tibble
#' @param log2fc_column a character defining the column name of the log2 foldchange values.
#' @param log2fc_min a numeric defining the minimum log2 foldchange to highlight.
#' @param significance_column a character defining the column name of the statistical significance values.
#' @param significance_max a numeric defining the maximum statistical significance to highlight.
#' @param labels_column a character defining the column name of the column for labeling.
#' @param show_lines a boolean for showing threshold lines.
#' @param point_size a numeric for changing the point size.
#' @param color_a a character defining the color for table_a expression.
#' @param color_b a character defining the color for table_b expression.
#' @param color_u a character defining the color for the union between both tables.
#'
#' @return a ggplot2 object
#' @export
#'
#' @examples
#' library(ggplot2, warn.conflicts = FALSE)
#' library(dplyr, warn.conflicts = FALSE)
#' library(tidyproteomics)
#'
#' # comparing two analytical methods, in substitute for two conditions
#' exp_a <- hela_proteins %>%
#'      expression(knockdown/control) %>%
#'      export_analysis(knockdown/control, .analysis = "expression")
#'
#' exp_b <- hela_proteins %>%
#'      expression(knockdown/control, .method = "limma") %>%
#'      export_analysis(knockdown/control, .analysis = "expression")
#'
#' plot_compexp(exp_a, exp_b, log2fc_min = 1, significance_column = "p_value") +
#'      ggplot2::labs(x = "(log2 FC) Wilcoxon Rank Sum",
#'                    y = "(log2 FC) Emperical Bayes (limma)",
#'                    title = "Hela p97 Knockdown ~ Control")
#'
plot_compexp <- function(
    table_a = NULL,
    table_b = NULL,
    log2fc_min = 2,
    log2fc_column = 'log2_foldchange',
    significance_max = 0.05,
    significance_column = 'adj_p_value',
    labels_column = 'protein',
    point_size = NULL,
    show_lines = TRUE,
    color_a = 'dodgerblue',
    color_b = 'firebrick1',
    color_u = 'purple'
){

  # visible bindings
  .data <- NULL

  tables <- compute_compexp(
    table_a,
    table_b,
    log2fc_min,
    log2fc_column,
    significance_max,
    significance_column,
    labels_column
  )

  table <- tables$merged
  table_grey <- tables$ignore
  table_a <- tables$sig_a_only
  table_b <- tables$sig_b_only
  table_i <- tables$sig_intersect
  table_label <- list(table_a, table_b, table_i) %>% dplyr::bind_rows()

  log2fc_column_a <- tables$log2fc_column_a
  log2fc_column_b <- tables$log2fc_column_b

  fc_min <- tables$fc_min
  fc_max <- tables$fc_max
  fc_scale <- tables$fc_scale

  signif_min <- tables$signif_min
  signif_max <- tables$signif_max

  # construct the volcano plot
  plot <- table %>%
    ggplot2::ggplot(ggplot2::aes(.data[[log2fc_column_a]], .data[[log2fc_column_b]]))

  # if(show_lines == TRUE){
  #   plot <- plot +
  #     ggplot2::geom_hline(yintercept = signif_max, color='black', linetype = 2, alpha=.33)

  if(fc_min != 0) {
    plot <- plot +
      ggplot2::geom_vline(xintercept = fc_min, color='black', linetype = 2, alpha=.33) +
      ggplot2::geom_vline(xintercept = -fc_min, color='black', linetype = 2, alpha=.33) +
      ggplot2::geom_hline(yintercept = fc_min, color='black', linetype = 2, alpha=.33) +
      ggplot2::geom_hline(yintercept = -fc_min, color='black', linetype = 2, alpha=.33)
  }
  # }

  plot <- plot +
    ggplot2::geom_point(data = table_grey, alpha=.5, color='grey') +
    ggplot2::geom_point(data = table_a, alpha=.5, color=color_a) +
    ggplot2::geom_point(data = table_b, alpha=.5, color=color_b) +
    ggplot2::geom_point(data = table_i, alpha=.5, color=color_u)

  # add the labels
  if(!is.null(labels_column)) {
    plot <- plot +
      ggrepel::geom_text_repel(data = table_label,
                               ggplot2::aes(label = .data[[labels_column]]),
                               hjust=-0.1, size=3)
  }

  # modify the color scheme
  plot <- plot +
    ggplot2::scale_y_continuous(breaks=fc_scale) +
    ggplot2::scale_x_continuous(breaks=fc_scale) +
    ggplot2::theme_bw() +
    # pretty up the axis labels
    ggplot2::labs(x = "Condition A",
                  y = "Condition B")

  return(plot)
}


#' Helper function to analysis between two expression tests
#'
#' @param table_a a tibble
#' @param table_b a tibble
#' @param log2fc_column a character defining the column name of the log2 foldchange values.
#' @param log2fc_min a numeric defining the minimum log2 foldchange to highlight.
#' @param significance_column a character defining the column name of the statistical significance values.
#' @param significance_max a numeric defining the maximum statistical significance to highlight.
#' @param labels_column a character defining the column name of the column for labeling.
#'
#' @return a list
#'
compute_compexp <- function(
    table_a = NULL,
    table_b = NULL,
    log2fc_min = 2,
    log2fc_column = 'log2_foldchange',
    significance_max = 0.05,
    significance_column = 'adj_p_value',
    labels_column = 'protein'
){

  # visible bindings
  .data <- NULL

  table <- list(
    table_a %>% dplyr::mutate(name = 'table_a'),
    table_b %>% dplyr::mutate(name = 'table_b')
  ) %>%
    dplyr::bind_rows() %>%
    dplyr::group_by_at(c('name', labels_column)) %>%
    dplyr::slice_min(.data[[significance_column]], n = 1, with_ties = F) %>%
    dplyr::ungroup()

  table_cols <- colnames(table)
  log2fc_column <- rlang::arg_match(log2fc_column, table_cols)
  significance_column <- rlang::arg_match(significance_column, table_cols)

  if(abs(mean(table[[log2fc_column]],na.rm=T) - 1) < 0.05) {
    cli::cli_div(theme = list(span.emph = list(color = "#ff4500"), span.info = list(color = "red")))
    cli::cli_alert_info(c("i" = "{.emph log2fc_column} ({.info {log2fc_column}}) appears to {.info NOT} be log transformed"))
  }

  table <- table %>%
    dplyr::filter(!is.na(.data[[log2fc_column]]),
                  !is.na(.data[[significance_column]])) %>%
    dplyr::mutate(keep = abs(.data[[log2fc_column]]) >= log2fc_min &
                    .data[[significance_column]] <= significance_max)

  if(length(which(table$keep == TRUE)) == 0){
    cli::cli_div(theme = list(span.emph = list(color = "#ff4500"), span.info = list(color = "red")))
    cli::cli_abort(c("i" = "no values significant at current settings",
                     "{.emph {log2fc_column}} ({.info {log2fc_min}}) and
                   {.emph {significance_column}} ({.info {significance_max}})"))
  }

  fc_min <- log2fc_min
  fc_max <- ceiling(max(abs(unlist(table[log2fc_column])), na.rm=T))
  fc_scale <- -5:5 * round(fc_max*2 / 11,1)

  signif_range <- table %>%
    dplyr::filter(.data[[significance_column]] <= significance_max) %>%
    dplyr::select(dplyr::all_of(significance_column)) %>%
    unlist() %>%
    range()

  signif_min <- min(signif_range) * .9
  signif_max <- max(signif_range) * 1.05

  table <- table %>%
    tidyr::pivot_wider(names_from = 'name', values_from = !dplyr::matches(labels_column))

  log2fc_column_a <- paste(log2fc_column, 'table_a', sep="_")
  log2fc_column_b <- paste(log2fc_column, 'table_b', sep="_")
  keep_column_a <- 'keep_table_a'
  keep_column_b <- 'keep_table_b'

  tables <- list(
    'merge' = table,
    'ignore' = table %>% dplyr::filter(.data[[keep_column_a]] == F) %>% dplyr::filter(.data[[keep_column_b]] == F),
    'sig_a_only' = table %>% dplyr::filter(.data[[keep_column_a]] == T) %>% dplyr::filter(.data[[keep_column_b]] == F),
    'sig_b_only' = table %>% dplyr::filter(.data[[keep_column_a]] == F) %>% dplyr::filter(.data[[keep_column_b]] == T),
    'sig_intersect' = table %>% dplyr::filter(.data[[keep_column_a]] == T) %>% dplyr::filter(.data[[keep_column_b]] == T),

    'log2fc_column_a' = log2fc_column_a,
    'log2fc_column_b' = log2fc_column_b,

    'fc_min' = fc_min,
    'fc_max' = fc_max,
    'fc_scale' = fc_scale,

    'signif_min' = signif_min,
    'signif_max' = signif_max
  )

  return(tables)
}

#' Comparative analysis between two expression tests
#'
#' @description
#' `export_compexp()` returns a table of the comparison in
#' expression differences between two methods or two sets of groups. For example,
#' one could run an expression difference for two different conditions (A and B)
#' prodived the experiment contained 3 samples condition A, condition B and WT,
#' then compare those results. The proteins showing up in the intersection
#' indicate common targets for condition A and B.
#' ```
#' expdiff_a <- protein_data %>%
#'    expression(experiment = "condition_a", control = "wt")
#'
#' expdiff_b <- protein_data %>%
#'    expression(experiment = "condition_b", control = "wt")
#'
#' export_compexp(expdiff_a, expdiff_b, export = "intersect")
#'
#' ```
#'
#' @param table_a a tibble
#' @param table_b a tibble
#' @param log2fc_column a character defining the column name of the log2 foldchange values.
#' @param log2fc_min a numeric defining the minimum log2 foldchange to highlight.
#' @param significance_column a character defining the column name of the statistical significance values.
#' @param significance_max a numeric defining the maximum statistical significance to highlight.
#' @param labels_column a character defining the column name of the column for labeling.
#' @param export a character string for the significance data to return
#'
#' @return a tibble
#' @export
#'
export_compexp <- function(
    table_a = NULL,
    table_b = NULL,
    log2fc_min = 2,
    log2fc_column = 'log2_foldchange',
    significance_max = 0.05,
    significance_column = 'adj_p_value',
    labels_column = 'protein',
    export = c('all', 'a_only', 'b_only', 'intersect')
){

  export <- rlang::arg_match(export)

  tables <- compute_compexp(
    table_a,
    table_b,
    log2fc_min,
    log2fc_column,
    significance_max,
    significance_column,
    labels_column
  )

  table_a <- tables$sig_a_only %>%
    dplyr::mutate(significant = 'a_only') %>%
    dplyr::select(dplyr::matches('protein|imputed|expression|foldchange|p_value|significant'))

  table_b <- tables$sig_b_only %>%
    dplyr::mutate(significant = 'b_only') %>%
    dplyr::select(dplyr::matches('protein|imputed|expression|foldchange|p_value|significant'))

  table_i <- tables$sig_intersect %>%
    dplyr::mutate(significant = 'ab_both') %>%
    dplyr::select(dplyr::matches('protein|imputed|expression|foldchange|p_value|significant'))

  table_label <- list(table_a, table_b, table_i) %>% dplyr::bind_rows()

  if(export == 'all') { return(table_label) }
  if(export == 'a_only') { return(table_a) }
  if(export == 'b_only') { return(table_b) }
  if(export == 'intersect') { return(table_i) }

  log2fc_column_a <- tables$log2fc_column_a
  log2fc_column_b <- tables$log2fc_column_b

  fc_min <- tables$fc_min
  fc_max <- tables$fc_max
  fc_scale <- tables$fc_scale

  signif_min <- tables$signif_min
  signif_max <- tables$signif_max

}
