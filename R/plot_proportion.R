#' Plot proportional expression values
#'
#' @description
#' `plot_proportion()` is a GGplot2 implementation for plotting the expression differences
#' as foldchange ~ scaled abundance. This allows for the visualization of selected
#' proteins See also `plot_volcano()`. This function can
#' take either a tidyproteomics data object or a table with the required headers.
#'
#' @param table a tibble
#' @param ... two sample comparison
#' @param log2fc_column a character defining the column name of the log2 foldchange values.
#' @param log2fc_min a numeric defining the minimum log2 foldchange to highlight.
#' @param significance_column a character defining the column name of the statistical significance values.
#' @param significance_max a numeric defining the maximum statistical significance to highlight.
#' @param proportion_column a character defining the column name of the proportional expression values.
#' @param proportion_min a numeric defining the minimum proportional expression to highlight.
#' @param labels_column a character defining the column name of the column for labeling.
#' @param label_significance a boolean for labeling values below the significance threshold.
#' @param show_pannels a boolean for showing colored up/down expression panels.
#' @param show_lines a boolean for showing threshold lines.
#' @param show_fc_scale a boolean for showing the secondary foldchange scale.
#' @param point_size a numeric for shanging the point size.
#' @param color_positive a character defining the color for positive (up) expression.
#' @param color_negative a character defining the color for negative (down) expression.
#' @param destination a character string
#' @param height a numeric
#' @param width a numeric

#'
#' @return a ggplot2 object
#' @export
#'
#' @examples
#' library(dplyr, warn.conflicts = FALSE)
#' library(tidyproteomics)
#' hela_proteins %>%
#'    expression(knockdown/control) %>%
#'    plot_proportion(knockdown/control, log2fc_min = 0.5, significance_column = 'p_value')
#'
#' # generates the same out come
#' # hela_proteins %>%
#' #    expression(knockdown/control) %>%
#' #    export_analysis(knockdown/control, .analysis = 'expression) %>%
#' #    plot_proportion(log2fc_min = 0.5, significance_column = 'p_value')
#'
#' # display the gene name instead
#' hela_proteins %>%
#'    expression(knockdown/control) %>%
#'    plot_proportion(knockdown/control, log2fc_min = 0.5, significance_column = 'p_value', labels_column = "gene_name")
#'
plot_proportion <- function(
    data = NULL,
    ...,
    log2fc_column = 'log2_foldchange',
    log2fc_min = 2,
    significance_column = 'adj_p_value',
    significance_max = 0.05,
    proportion_column = 'proportional_expression',
    proportion_min = 0.01,
    labels_column = 'protein',
    label_significance = TRUE,
    show_pannels = FALSE,
    show_lines = TRUE,
    show_fc_scale = TRUE,
    point_size = NULL,
    color_positive = 'dodgerblue',
    color_negative = 'firebrick1',
    destination = 'plot',
    height = 5,
    width = 8
){

  # visible bindings
  .data <- NULL
  metric <- NULL

  file_name <- 'proportion'
  if('tidyproteomics' %in% class(data)) {
    table <- data %>% export_analysis(..., .analysis = 'expression')
    str_quo <- tidyproteomics_quo_name(...)
    file_name = glue::glue("{data$analyte}_{data$quantitative_source}_proportion_{str_quo}")
  } else {
    table <- data
  }

  table_cols <- colnames(table)
  log2fc_column <- rlang::arg_match(log2fc_column, table_cols)
  significance_column <- rlang::arg_match(significance_column, table_cols)

  if(abs(mean(table[[log2fc_column]],na.rm=T) - 1) < 0.05) {
    cli::cli_div(theme = list(span.emph = list(color = "#ff4500"), span.info = list(color = "red")))
    cli::cli_alert_info(c("i" = "{.emph log2fc_column} ({.info {log2fc_column}}) appears to {.info NOT} be log transformed"))
  }

  if(abs(sum(table[[proportion_column]], na.rm=T) - 1) <= 0.05){
    cli::cli_div(theme = list(span.emph = list(color = "#ff4500"), span.info = list(color = "blue")))
    cli::cli_alert_info(c("{.info {proportion_column}} appears to sum to 1",
                          " adjusting values to 100(%)"))

    proportion_min <- proportion_min * 100
    table[[proportion_column]] = table[[proportion_column]] * 100

  } else if(abs(sum(table[[proportion_column]], na.rm=T) - 100) <= 5){
    cli::cli_div(theme = list(span.emph = list(color = "#ff4500"), span.info = list(color = "blue")))
    cli::cli_alert_info(c("{.info {proportion_column}} appears to sum to 100, assuming (%)"))
  } else {
    cli::cli_div(theme = list(span.emph = list(color = "#ff4500"), span.info = list(color = "blue")))
    cli::cli_alert_info(c("{.info {proportion_column}} does not appear to sum to 1 or 100"))
  }

  table <- table %>%
    dplyr::filter(!is.na(.data[[log2fc_column]]),
                  !is.na(.data[[significance_column]])) %>%
    dplyr::mutate(keep = abs(.data[[log2fc_column]]) >= log2fc_min &
                    .data[[significance_column]] <= significance_max,
                  show = .data[[proportion_column]] >= proportion_min)

  if(length(which(table$keep == TRUE)) == 0){
    show_pannels = FALSE
    show_lines = FALSE

    cli::cli_div(theme = list(span.emph = list(color = "#ff4500"), span.info = list(color = "red")))
    cli::cli_inform(c("i" = "no values significant at current settings",
                      "i" = "... {.emph {log2fc_column}} ({.info {log2fc_min}}) and
                   {.emph {significance_column}} ({.info {significance_max}})"))

    if(min(table[,significance_column]) == 1) {
      significance_column <- 'p_value'
      cli::cli_inform(c("i" = "... using {.emph {significance_column}}"))
    }

    table <- table %>%
      dplyr::mutate(metric = abs(.data[[log2fc_column]]) + 10 * (1 - .data[[significance_column]])) %>%
      dplyr::arrange(dplyr::desc(metric)) %>%
      dplyr::mutate(rank = dplyr::row_number()) %>%
      dplyr::mutate(keep = rank <= 10)
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
  proportion_max <- max(table[[proportion_column]])

  table_show <- table %>% dplyr::filter(.data[['show']] == T)
  table_grey <- table %>% dplyr::filter(.data[['keep']] == F)
  table_label <- table %>% dplyr::filter(.data[['keep']] == T)
  table_neg <- table_label %>% dplyr::filter(.data[[log2fc_column]] < 0)
  table_pos <- table_label %>% dplyr::filter(.data[[log2fc_column]] > 0)

  # construct the volcano plot
  if(is.null(point_size)) {
    plot <- table %>%
      ggplot2::ggplot(ggplot2::aes(.data[[log2fc_column]], .data[[proportion_column]]))
  } else {
    point_size <- rlang::arg_mathc(point_size, table_cols)
    plot <- table %>%
      ggplot2::ggplot(ggplot2::aes(.data[[log2fc_column]], .data[[proportion_column]],
                                   size = .data[[point_size]]))
  }

  if(show_lines == TRUE){
    plot <- plot +
      ggplot2::geom_hline(yintercept = proportion_min, color='black', linetype = 2, alpha=.33)

    if(fc_min != 0) {
      plot <- plot +
        ggplot2::geom_vline(xintercept = fc_min, color='black', linetype = 2, alpha=.33) +
        ggplot2::geom_vline(xintercept = -fc_min, color='black', linetype = 2, alpha=.33)
    }
  }

  plot <- plot +
    ggplot2::geom_point(data = table_grey, alpha=.5, color='grey') +
    ggplot2::geom_point(data = table_neg, alpha=.5, color=color_negative) +
    ggplot2::geom_point(data = table_pos, alpha=.5, color=color_positive)

  if(show_pannels == TRUE) {
    plot <- plot +
      ggplot2::annotate("rect",
                        xmin = -Inf, xmax = -fc_min,
                        ymin = proportion_min, ymax = proportion_max,
                        alpha = .1, fill=color_negative) +

      ggplot2::annotate("rect",
                        xmin = fc_min, xmax = Inf,
                        ymin = proportion_min, ymax = proportion_max,
                        alpha = .1, fill=color_positive)
  }

  # add the labels
  if(!is.null(labels_column)) {
    plot <- plot +
      ggrepel::geom_text_repel(data = table_show,
                               ggplot2::aes(label = .data[[labels_column]]),
                               hjust=-0.1, size=3)

    if(label_significance == TRUE) {
      plot <- plot +
        ggrepel::geom_text_repel(data = table_label,
                                 ggplot2::aes(label = .data[[labels_column]]),
                                 hjust=-0.1, size=3)
    }
  }

  # modify the color scheme
  plot <- plot +
    ggplot2::scale_color_manual(values = theme_palette()) +
    ggplot2::scale_y_log10(breaks=c(30,10,5,3,1,c(1,3,5)/10, c(1,5)/100)) +
    ggplot2::scale_x_continuous(breaks=fc_scale) +
    ggplot2::theme_bw() +
    # pretty up the axis labels
    ggplot2::labs(x = paste0("Expression Difference (",log2fc_column, ")"),
                  y = paste0("Proportional Expression (%)"))

  if(show_fc_scale == TRUE) {
    plot <- plot +
      ggplot2::annotate('text',
                        y=Inf, x=fc_scale, size=3,
                        label = round(2^abs(fc_scale),1), vjust=-1) +
      ggplot2::annotate('text',
                        y=Inf, x=-Inf, size=3,
                        label = 'FoldChange', vjust=-2.5, hjust=-.5)
  }

  return(plot_save(plot,
                   data,
                   file_name = file_name,
                   destination = destination,
                   height = height,
                   width = width))
}
