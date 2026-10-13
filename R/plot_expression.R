#' Plot protein expression differences
#'
#' @description
#' `plot_expression()` is a ggplot2 implementation for visualizing protein abundance
#' distributions across two sample conditions for top differentially expressed proteins.
#' The resulting plot is faceted by protein/annotation label, with the two samples
#' on the x-axis, log2 abundance on the y-axis, a boxplot representing the distribution,
#' and individual replicate points overlaid (filled circle for non-imputed, hollow circle
#' for imputed).
#'
#' @param data a tidyproteomics data object
#' @param ... two sample comparison expression (e.g. `knockdown/control`)
#' @param log2fc_min a numeric defining the minimum log2 foldchange threshold. Default is 1.
#' @param log2fc_column a character defining the column name of log2 foldchange values.
#' @param significance_max a numeric defining the maximum significance threshold. Default is 0.05.
#' @param significance_column a character defining the column name of significance values.
#' @param top_n an integer defining the top N proteins to plot post-filtering. Default is 10.
#' @param regulation a character defining direction: `"both"` (default), `"up"`, or `"down"`.
#' @param labels_column a character defining the annotation column for facet titles (default is `"protein"`).
#' @param show_title input FALSE, TRUE for an auto-generated title or any character string.
#' @param destination a character string (`"plot"`, `"save"`, `"png"`, `"svg"`, `"tiff"`, `"jpeg"`).
#' @param height a numeric plot height in inches. Default is 6.
#' @param width a numeric plot width in inches. Default is 8.
#'
#' @return a ggplot2 object or tidyproteomics data object depending on `destination`
#' @export
#'
#' @examples
#' library(dplyr, warn.conflicts = FALSE)
#' library(tidyproteomics)
#'
#' hela_proteins %>%
#'   expression(knockdown/control) %>%
#'   plot_expression(knockdown/control, log2fc_min = 0.5, significance_column = "p_value", top_n = 6)
#'
plot_expression <- function(
    data = NULL,
    ...,
    log2fc_min = 1,
    log2fc_column = 'log2_foldchange',
    significance_max = 0.05,
    significance_column = 'adj_p_value',
    top_n = 10,
    regulation = c('both', 'up', 'down'),
    labels_column = 'protein',
    show_title = TRUE,
    destination = 'plot',
    height = 6,
    width = 8
){

  # visible bindings
  .data <- NULL
  imputed <- NULL
  imputed_status <- NULL
  sample <- NULL
  plot_label <- NULL

  check_data(data)
  regulation <- match.arg(regulation)

  str_quo <- tidyproteomics_quo(...)
  if(is.null(str_quo)) {
    cli::cli_alert_warning("No comparison provided")
    cli::cli_alert_info("Try one of ({paste(names(data$analysis), collapse = ', ')})")
    return(data)
  }

  sample_a <- str_quo[['variable']]
  sample_b <- str_quo[['value']]
  if(str_quo['operator'] != "/") {
    cli::cli_div(theme = list(span.emph = list(color = "#ff4500")))
    cli::cli_abort("Comparison operator must be {.emph \"/\"} (e.g. {sample_a}{.emph /}{sample_b})")
  }
  check_samples(data, sample_a, sample_b)

  str_quo_name <- tidyproteomics_quo_name(...)
  file_name <- glue::glue("{data$analyte}_{data$quantitative_source}_expression_{str_quo_name}")

  table <- data %>% export_analysis(..., .analysis = 'expression')
  table_cols <- colnames(table)

  log2fc_column <- rlang::arg_match(log2fc_column, table_cols)
  significance_column <- rlang::arg_match(significance_column, table_cols)

  available_labels <- unique(c(data$identifier, get_annotation_terms(data), table_cols))
  labels_column <- rlang::arg_match(labels_column, available_labels)

  # filter by foldchange, significance, and regulation direction
  table <- table %>%
    dplyr::filter(!is.na(.data[[log2fc_column]]),
                  !is.na(.data[[significance_column]]),
                  !is.infinite(.data[[log2fc_column]]))

  table_filtered <- switch(
    regulation,
    'up' = table %>%
      dplyr::filter(.data[[log2fc_column]] >= log2fc_min,
                    .data[[significance_column]] <= significance_max) %>%
      dplyr::arrange(dplyr::desc(.data[[log2fc_column]])),
    'down' = table %>%
      dplyr::filter(.data[[log2fc_column]] <= -log2fc_min,
                    .data[[significance_column]] <= significance_max) %>%
      dplyr::arrange(.data[[log2fc_column]]),
    'both' = table %>%
      dplyr::filter(abs(.data[[log2fc_column]]) >= log2fc_min,
                    .data[[significance_column]] <= significance_max) %>%
      dplyr::arrange(dplyr::desc(abs(.data[[log2fc_column]])))
  )

  if(nrow(table_filtered) == 0){
    cli::cli_div(theme = list(span.emph = list(color = "#ff4500"), span.info = list(color = "red")))
    cli::cli_inform(c("i" = "no values significant at current settings",
                      "i" = "... {.emph {log2fc_column}} ({.info {log2fc_min}}) and {.emph {significance_column}} ({.info {significance_max}})",
                      "i" = "... falling back to top candidates ranked by effect and significance"))

    table_filtered <- table %>%
      dplyr::mutate(metric = abs(.data[[log2fc_column]]) + 10 * (1 - .data[[significance_column]])) %>%
      dplyr::arrange(dplyr::desc(metric))
  }

  table_top <- table_filtered %>% head(top_n)

  # establish facet labels preserving rank order
  raw_labels <- table_top[[labels_column]]
  labels_clean <- ifelse(is.na(raw_labels) | trimws(as.character(raw_labels)) == "",
                         table_top[[data$identifier[1]]],
                         as.character(raw_labels))
  labels_clean <- stringr::str_extract(labels_clean, "^[^;]+(;\\s*[^;]+)?")
  labels_unique <- make.unique(labels_clean)

  table_top$plot_label <- factor(labels_unique, levels = labels_unique)
  mapping_table <- table_top %>%
    dplyr::select(dplyr::all_of(c(data$identifier[1], "plot_label"))) %>%
    dplyr::distinct()

  # extract long-format quantitative data
  quant_source <- data$quantitative_source
  quant_col <- paste0("abundance_", quant_source)

  long_data <- as.data.frame(data, shape = 'long') %>%
    dplyr::filter(.data[[data$identifier[1]]] %in% table_top[[data$identifier[1]]],
                  sample %in% c(sample_a, sample_b),
                  !is.na(.data[[quant_col]])) %>%
    dplyr::left_join(mapping_table, by = data$identifier[1])

  # encode imputed status
  long_data <- long_data %>%
    dplyr::mutate(
      imputed_status = factor(
        ifelse(as.character(imputed) %in% c("1", "TRUE"), "imputed", "observed"),
        levels = c("observed", "imputed")
      ),
      sample = factor(sample, levels = c(sample_a, sample_b)),
      log2_abundance = log2(.data[[quant_col]])
    )

  # calculate facet-specific annotation heights based on data span
  table_annot <- long_data %>%
    dplyr::group_by(plot_label) %>%
    dplyr::summarise(
      min_y = min(log2_abundance, na.rm = TRUE),
      max_y = max(log2_abundance, na.rm = TRUE),
      .groups = 'drop'
    ) %>%
    dplyr::inner_join(table_top, by = "plot_label") %>%
    dplyr::mutate(
      y_span = pmax(max_y - min_y, 0.5),
      y_bar = max_y + y_span * 0.08,
      y_text = max_y + y_span * 0.16,
      is_sig = ifelse(p_value < 0.05, "*", ""),
      is_sig = ifelse(p_value < 0.01, "**", is_sig)
    )

  print(table_annot)

  # build title
  title_txt <- str_quo_name %>%
    stringr::str_replace("-", " / ") %>%
    stringr::str_replace_all("_", " ") %>%
    stringr::str_to_title()
  if(is.character(show_title)) {
    title_txt <- show_title
  }

  subtitle_txt <- glue::glue("Top {nrow(table_top)} ({regulation}) | {significance_column} <= {significance_max}, |{log2fc_column}| >= {log2fc_min}")
  if(isFALSE(show_title)) {
    title_txt <- NULL
    subtitle_txt <- NULL
  }

  # generate plot
  plot <- ggplot2::ggplot(long_data, ggplot2::aes(x = sample, y = log2_abundance)) +
    ggplot2::geom_boxplot(ggplot2::aes(fill = sample),
                          color = 'grey40',
                          outlier.shape = NA,
                          alpha = 0.5) +
    ggplot2::geom_point(
      ggplot2::aes(shape = imputed_status),
      position = ggplot2::position_jitter(width = 0.15, height = 0),
      size = 2.5,
      alpha = 0.85
    ) +
    ggplot2::scale_shape_manual(
      values = c("observed" = 16, "imputed" = 1),
      drop = FALSE,
      name = "Status"
    ) +
    ggplot2::scale_fill_manual(
      values = theme_palette(2),
      name = "Sample"
    ) +
    ggplot2::facet_wrap(~ plot_label, scales = "free_y") +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.05, 0.22))) +
    ggplot2::theme_classic() +
    ggplot2::labs(
      x = "Sample",
      y = expression(Log[2]~Abundance),
      title = title_txt,
      subtitle = subtitle_txt
    ) +
    ggplot2::theme(
      legend.position = "none"
    ) +
    ggplot2::geom_segment(
      data = table_annot,
      ggplot2::aes(x = 1, xend = 2, y = y_bar, yend = y_bar),
      color = 'red',
      inherit.aes = FALSE
    ) +
    ggplot2::geom_text(
      data = table_annot,
      ggplot2::aes(x = 1.5, y = y_text, label = paste(is_sig, "p =", signif(.data[[significance_column]], 2))),
      color = 'red',
      vjust = 0,
      hjust = 0.5,
      size = 3,
      inherit.aes = FALSE
    )

  return(plot_save(
    plot,
    data,
    file_name = file_name,
    destination = destination,
    height = height,
    width = width
  ))
}
