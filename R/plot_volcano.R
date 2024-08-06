#' Volcano plot of expression values
#'
#' @description
#' `plot_volcano()` is a GGplot2 implementation for plotting the expression differences
#' as foldchange ~ statistical significance. See also `plot_proportion()`.  This function can
#' take either a tidyproteomics data object or a table with the required headers.
#'
#' @param data a tibble
#' @param ... two sample comparison
#' @param log2fc_column a character defining the column name of the log2 foldchange values.
#' @param log2fc_min a numeric defining the minimum log2 foldchange to highlight.
#' @param significance_column a character defining the column name of the statistical significance values.
#' @param significance_max a numeric defining the maximum statistical significance to highlight.
#' @param labels_column a character defining the column name of the column for labeling.
#' @param show_pannels a boolean for showing colored up/down expression panels.
#' @param show_lines a boolean for showing threshold lines.
#' @param show_fc_scale a boolean for showing the secondary foldchange scale.
#' @param show_title input FALSE, TRUE for an auto-generated title or any charcter string.
#' @param show_pval_1 a boolean for showing expressions with pvalue == 1.
#' @param point_size a character reference to a numerical value in the expression table
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
#'    plot_volcano(knockdown/control, log2fc_min = 0.5, significance_column = "p_value")
#'
#' # generates the same out come
#' # hela_proteins %>%
#' #     expression(knockdown/control) %>%
#' #     export_analysis(knockdown/control, .analysis = "expression") %>%
#' #     plot_volcano(log2fc_min = 0.5, significance_column = "p_value")
#'
#' # display the gene name instead
#' hela_proteins %>%
#'    expression(knockdown/control) %>%
#'    plot_volcano(knockdown/control, log2fc_min = 0.5, significance_column = "p_value", labels_column = "gene_name")
#'
plot_volcano <- function(
    data = NULL,
    ...,
    log2fc_min = 1,
    log2fc_column = 'log2_foldchange',
    significance_max = 0.05,
    significance_column = 'adj_p_value',
    labels_column = 'gene_name',
    show_pannels = TRUE,
    show_lines = TRUE,
    show_fc_scale = TRUE,
    show_title = TRUE,
    show_pval_1 = TRUE,
    point_size = NULL,
    color_positive = 'dodgerblue',
    color_negative = 'firebrick1',
    destination = 'plot',
    height = 5,
    width = 8
){

  # visible bindings
  metric <- NULL
  significance_plot <- 'p_value'
  title_txt <- NULL

  file_name = "volcano"
  if('tidyproteomics' %in% class(data)) {
    table <- data %>% export_analysis(..., .analysis = 'expression')
    str_quo <- tidyproteomics_quo_name(...)
    title_txt <- str_quo %>%
      stringr::str_replace("-", " / ") %>%
      stringr::str_replace_all("_", " ") %>%
      stringr::str_to_title()


    file_name = glue::glue("{data$analyte}_{data$quantitative_source}_volcano_{str_quo}")
  } else {
    table <- data
  }

  if(is.character(show_title)){title_txt <- show_title}

  table_cols <- colnames(table)
  log2fc_column <- rlang::arg_match(log2fc_column, table_cols)
  significance_column <- rlang::arg_match(significance_column, table_cols)
  if(!is.null(labels_column)) {
    labels_column <- rlang::arg_match(labels_column, table_cols)
  }

  if(abs(mean(table[[log2fc_column]],na.rm=T) - 1) < 0.05) {
    cli::cli_div(theme = list(span.emph = list(color = "#ff4500"), span.info = list(color = "red")))
    cli::cli_alert_info(c("i" = "{.emph log2fc_column} ({.info {log2fc_column}}) appears to {.info NOT} be log transformed"))
  }

  table <- table %>%
    dplyr::filter(!is.na(.data[[log2fc_column]])) %>%
    dplyr::filter(!is.na(.data[[significance_column]])) %>%
    dplyr::filter(!is.infinite(.data[[log2fc_column]])) %>%
    dplyr::mutate(keep = abs(.data[[log2fc_column]]) >= log2fc_min &
                    .data[[significance_column]] <= significance_max)

  title_n_exp <- table %>% nrow()

  if(show_pval_1 == FALSE){
    table <- table %>% dplyr::filter(!.data[['p_value']] == 1)
  }

  show_signif <- TRUE
  if(length(which(table$keep == TRUE)) == 0){
    show_pannels = FALSE
    show_lines = FALSE
    show_signif <- FALSE

    cli::cli_div(theme = list(span.emph = list(color = "#ff4500"), span.info = list(color = "red")))
    cli::cli_inform(c("i" = "no values significant at current settings",
                      "i" = "... {.emph {log2fc_column}} ({.info {log2fc_min}}) and
                   {.emph {significance_column}} ({.info {significance_max}})"))

    if(min(table[,significance_column]) == 1) {
      significance_column <- 'p_value'
      cli::cli_inform(c("i" = "... using {.emph {significance_column}} as significance value"))
    }

    table <- table %>%
      dplyr::mutate(metric = abs(.data[[log2fc_column]]) + 10 * (1 - .data[[significance_column]])) %>%
      dplyr::arrange(dplyr::desc(metric)) %>%
      dplyr::mutate(rank = dplyr::row_number()) %>%
      dplyr::mutate(keep = rank <= 10)
  }

  fc_min <- log2fc_min
  fc_max <- ceiling(max(abs(unlist(table[,log2fc_column])), na.rm=T))
  fc_scale <- -5:5 * round(fc_max*2 / 11,1)

  signif_range <- table %>%
    dplyr::filter(.data[[significance_column]] <= significance_max) %>%
    dplyr::select(dplyr::all_of(significance_plot)) %>%
    unlist() %>%
    range()

  signif_min <- min(signif_range) * .9
  signif_max <- max(signif_range) * 1.05

  table_grey <- table %>% dplyr::filter(.data[['keep']] == F)
  table_label <- table %>% dplyr::filter(.data[['keep']] == T)
  table_neg <- table_label %>% dplyr::filter(.data[[log2fc_column]] < 0)
  table_pos <- table_label %>% dplyr::filter(.data[[log2fc_column]] > 0)

  # construct the volcano plot
  if(is.null(point_size)) {
    plot <- table %>%
      ggplot2::ggplot(ggplot2::aes(.data[[log2fc_column]], .data[[significance_plot]]))
  } else {
    point_size <- rlang::arg_match(point_size, table_cols)
    plot <- table %>%
      ggplot2::ggplot(ggplot2::aes(.data[[log2fc_column]], .data[[significance_plot]],
                                   size = .data[[point_size]]))
  }

  if(show_lines == TRUE){
    plot <- plot +
      ggplot2::geom_hline(yintercept = signif_max, color='black', linetype = 2, alpha=.33)

    if(fc_min != 0) {
      plot <- plot +
        ggplot2::geom_vline(xintercept = fc_min, color='black', linetype = 2, alpha=.33) +
        ggplot2::geom_vline(xintercept = -fc_min, color='black', linetype = 2, alpha=.33)
    }
  }

  if(show_signif == TRUE){
    plot <- plot +
      ggplot2::geom_point(data = table_grey, alpha=.5, color='grey') +
      ggplot2::geom_point(data = table_neg, alpha=.5, color=color_negative) +
      ggplot2::geom_point(data = table_pos, alpha=.5, color=color_positive)
  } else {
    plot <- plot +
      ggplot2::geom_point(data = table, alpha=.5, color='grey') +
      ggplot2::geom_point(data = table_neg, shape = 1) +
      ggplot2::geom_point(data = table_pos, shape = 1)
  }

  if(show_pannels == TRUE) {
    plot <- plot +
      ggplot2::annotate("rect",
                        xmin = -Inf, xmax = -fc_min,
                        ymin = signif_min, ymax = signif_max,
                        alpha = .1, fill=color_negative) +

      ggplot2::annotate("rect",
                        xmin = fc_min, xmax = Inf,
                        ymin = signif_min, ymax = signif_max,
                        alpha = .1, fill=color_positive)
  }

  # add the labels
  if(!is.null(labels_column)) {
    plot <- plot +
      ggrepel::geom_text_repel(data = table_label,
                               ggplot2::aes(label = .data[[labels_column]]),
                               hjust=-0.1, size=3)
  }

  significance_plot_txt <- significance_plot %>% stringr::str_replace_all("_"," ") %>% stringr::str_to_title()
  log2fc_column_txt <- log2fc_column %>% stringr::str_replace_all("_"," ") %>% stringr::str_to_title()

  n_sig <- nrow(table_label)
  if(show_signif == FALSE){ n_sig <- 0 }
  subtitle_txt <- glue::glue("{n_sig} of {title_n_exp} found significant: {significance_column} <= {significance_max}, {log2fc_column} >= {log2fc_min}")
  if(show_title == FALSE){ subtitle_txt <- NULL; title_txt <- NULL }

  # modify the color scheme
  plot <- plot +
    ggplot2::scale_color_manual(values = theme_palette()) +
    ggplot2::scale_y_continuous(trans=reverselog_transformation(10), breaks = signif(1/10^(0:20), 1)) +
    ggplot2::scale_x_continuous(breaks=fc_scale) +
    ggplot2::geom_hline(yintercept = 1, color = NA) +
    ggplot2::theme_classic() +
    # pretty up the axis labels
    ggplot2::labs(x = glue::glue("Expression Difference ({log2fc_column_txt})"),
                  y = glue::glue("Significance ({significance_plot_txt})"),
                  title = title_txt,
                  subtitle = subtitle_txt) +
    ggplot2::theme(
      plot.subtitle = ggplot2::element_text(size=9, face="italic"),
      axis.title = ggplot2::element_text(size=10)
    )

  if(show_fc_scale == TRUE) {
    plot <- plot +
      ggplot2::annotate('text',
                        y=Inf, x=fc_scale, size=3, color='grey40',
                        label = round(2^abs(fc_scale),1), vjust=-1) +
      ggplot2::annotate('text',
                        y=Inf, x=-Inf, size=3, , color='grey40',
                        label = 'FoldChange', vjust=-2.5, hjust=-.5)
  }

  return(plot_save(plot,
                   data,
                   file_name = file_name,
                   destination = destination,
                   height = height,
                   width = width))
}
