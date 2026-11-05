#' Plot Candidate Targets and Effect Sizes
#'
#' @description
#' `plot_candidates()` visualizes effect sizes (log2 foldchange) along with Wald-type
#' confidence intervals for specific candidate proteins or literature-nominated
#' validation targets. Candidates can optionally be categorized and colored by validation status.
#'
#' @param data a tidyproteomics data object (or an exported expression table)
#' @param ... two sample comparison expression (e.g. `knockdown/control`)
#' @param candidates a character vector of identifiers or gene names to include in the plot.
#'   If NULL, defaults to `proteins` or the top 15 candidates ranked by statistical significance and effect size.
#' @param proteins alias for `candidates` for backward compatibility.
#' @param labels_column a character defining the column name for feature labeling. Default is "gene_name".
#' @param status an optional named character vector of validation status tags corresponding to gene/protein
#'   identifiers (e.g. `c(HERPUD1 = "Validated", RACGAP1 = "Empirical change")`).
#' @param status_column an optional column name in the expression table containing status annotations.
#' @param confidence numeric confidence level for intervals. Default is 0.95 (95% CI).
#' @param log2fc_column a character defining the column name of log2 foldchange values. Default is "log2_foldchange".
#' @param significance_column a character defining the column name of significance values. Default is "adj_p_value".
#' @param significance_max a numeric defining the significance threshold. Default is 0.05.
#' @param show_title boolean or character string for plot title. Default is TRUE.
#' @param tag optional character string to append to the plot file name.
#' @param destination character string indicating output type: "plot", "save", "png", "svg", "tiff", "jpeg". Default is "plot".
#' @param height numeric plot height in inches. Default is 6.
#' @param width numeric plot width in inches. Default is 8.
#' @param path optional path to save destination plot.
#' @param ... passthrough ggplot2::ggsave arguments.
#'
#' @return a ggplot2 object or tidyproteomics data object depending on `destination`.
#' @export
#'
#' @examples
#' library(dplyr, warn.conflicts = FALSE)
#' library(tidyproteomics)
#'
#' hela_proteins %>%
#'   expression(knockdown/control) %>%
#'   plot_candidates(
#'     knockdown/control,
#'     candidates = c("RPL26", "RPS27A", "CKS2"),
#'     labels_column = "gene_name"
#'   )
#'
plot_candidates <- function(
    data = NULL,
    ...,
    candidates = NULL,
    proteins = NULL,
    labels_column = "gene_name",
    status = NULL,
    status_column = NULL,
    confidence = 0.95,
    log2fc_column = "log2_foldchange",
    significance_column = "adj_p_value",
    significance_max = 0.05,
    show_title = TRUE,
    tag = NULL,
    destination = "plot",
    height = 6,
    width = 8,
    path = NULL
){

  # visible bindings
  .data <- NULL
  plot_label <- NULL
  ci_lower <- NULL
  ci_upper <- NULL

  # extract comparison and expression data
  if ("tidyproteomics" %in% class(data)) {
    check_data(data)
    str_quo <- tidyproteomics_quo(...)
    if (is.null(str_quo)) {
      cli::cli_alert_warning("No comparison provided")
      cli::cli_alert_info("Try one of ({paste(names(data$analysis), collapse = ', ')})")
      return(data)
    }
    sample_a <- str_quo[['variable']]
    sample_b <- str_quo[['value']]
    if (str_quo['operator'] != "/") {
      cli::cli_div(theme = list(span.emph = list(color = "#ff4500")))
      cli::cli_abort("Comparison operator must be {.emph \"/\"} (e.g. {sample_a}{.emph /}{sample_b})")
    }
    check_samples(data, sample_a, sample_b)
    str_quo_name <- tidyproteomics_quo_name(...)
    table <- data %>% export_analysis(..., .analysis = 'expression')
    analyte <- data$analyte
    quant_source <- data$quantitative_source
    id_col <- data$identifier[1]
  } else if (is.data.frame(data)) {
    table <- data
    str_quo_name <- "expression"
    analyte <- "proteins"
    quant_source <- "normalized"
    id_col <- if ("protein" %in% colnames(table)) "protein" else colnames(table)[1]
  } else {
    cli::cli_abort("`data` must be a tidyproteomics object or data.frame")
  }

  table_cols <- colnames(table)
  log2fc_column <- rlang::arg_match(log2fc_column, table_cols)
  significance_column <- rlang::arg_match(significance_column, table_cols)

  available_labels <- unique(c(id_col, table_cols))
  if (labels_column %in% available_labels) {
    label_var <- labels_column
  } else {
    label_var <- id_col
  }

  # resolve candidates alias with proteins
  if (is.null(candidates) && !is.null(proteins)) {
    candidates <- proteins
  }

  # filter to nominated candidates or fallback to top candidates
  if (!is.null(candidates) && length(candidates) > 0) {
    targets <- as.character(candidates)
    sub <- table %>%
      dplyr::filter(
        (!is.na(.data[[label_var]]) & .data[[label_var]] %in% targets) |
          (!is.na(.data[[id_col]]) & .data[[id_col]] %in% targets)
      )
  } else {
    # fallback to top 15 candidates ranked by significance and foldchange
    sub <- table %>%
      dplyr::filter(!is.na(.data[[log2fc_column]]), !is.na(.data[[significance_column]])) %>%
      dplyr::mutate(rank_metric = abs(.data[[log2fc_column]]) + 5 * (1 - pmin(.data[[significance_column]], 1))) %>%
      dplyr::arrange(dplyr::desc(rank_metric)) %>%
      head(15)
  }

  if (nrow(sub) == 0) {
    cli::cli_alert_warning("No matching candidates found in data.")
    return(if ("tidyproteomics" %in% class(data)) data else NULL)
  }

  # clean labels
  raw_labels <- sub[[label_var]]
  labels_clean <- ifelse(is.na(raw_labels) | trimws(as.character(raw_labels)) == "",
                         sub[[id_col]],
                         as.character(raw_labels))
  labels_clean <- stringr::str_extract(labels_clean, "^[^;]+(;\\s*[^;]+)?")
  sub$plot_label <- labels_clean

  # compute standard error and confidence intervals
  fc <- as.numeric(sub[[log2fc_column]])

  if ("limma_t_statistic" %in% colnames(sub) && any(!is.na(sub$limma_t_statistic))) {
    t_stat <- as.numeric(sub$limma_t_statistic)
    se <- ifelse(!is.na(t_stat) & abs(t_stat) > 1e-4, abs(fc / t_stat), NA_real_)
  } else {
    se <- rep(NA_real_, nrow(sub))
  }

  # fallback approximation from p-value if se is not available
  p_col <- if ("p_value" %in% colnames(sub)) "p_value" else significance_column
  p_vals <- pmin(pmax(as.numeric(sub[[p_col]]), 1e-15), 0.95)
  z_vals <- pmax(stats::qnorm(1 - p_vals / 2), 0.05)
  se_from_p <- abs(fc / z_vals)
  se <- ifelse(is.na(se) | is.infinite(se) | se <= 0, se_from_p, se)
  se <- pmin(se, 10)

  alpha <- 1 - confidence
  z_crit <- stats::qnorm(1 - alpha / 2)
  sub$ci_lower <- fc - z_crit * se
  sub$ci_upper <- fc + z_crit * se

  # resolve status annotation
  if (!is.null(status) && length(status) > 0) {
    sub$status <- ifelse(
      sub[[label_var]] %in% names(status),
      status[sub[[label_var]]],
      ifelse(sub[[id_col]] %in% names(status), status[sub[[id_col]]], "Nominated")
    )
  } else if (!is.null(status_column) && status_column %in% colnames(sub)) {
    sub$status <- as.character(sub[[status_column]])
  } else {
    sub$status <- ifelse(sub[[significance_column]] <= significance_max, "Significant", "Non-significant")
  }

  # order rows by log2 foldchange ascending (so highest foldchange is at top of plot)
  sub <- sub %>% dplyr::arrange(.data[[log2fc_column]])
  sub$plot_label <- factor(sub$plot_label, levels = unique(sub$plot_label))

  # color palette for status
  statuses <- unique(sub$status)
  num_statuses <- length(statuses)
  status_palette <- theme_palette(max(num_statuses, 2))[1:num_statuses]
  names(status_palette) <- statuses

  # title and subtitle
  title_txt <- glue::glue("{stringr::str_to_title(gsub('_', ' ', str_quo_name))} — Candidate Targets")
  if (is.character(show_title)) {
    title_txt <- show_title
  } else if (isFALSE(show_title)) {
    title_txt <- NULL
  }

  conf_pct <- round(confidence * 100)
  subtitle_txt <- glue::glue("Effect sizes with {conf_pct}% Wald CI | {nrow(sub)} candidates")
  if (isFALSE(show_title)) {
    subtitle_txt <- NULL
  }

  # construct plot
  plot <- ggplot2::ggplot(
    sub,
    ggplot2::aes(x = .data[[log2fc_column]], y = plot_label, color = status)
  ) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "grey50", alpha = 0.7) +
    ggplot2::geom_errorbar(
      ggplot2::aes(xmin = ci_lower, xmax = ci_upper),
      width = 0.25,
      linewidth = 0.7
    ) +
    ggplot2::geom_point(size = 3) +
    ggplot2::scale_color_manual(values = status_palette, name = "Status") +
    ggplot2::theme_classic() +
    ggplot2::theme(
      axis.text.y = ggplot2::element_text(size = 9, color = "grey20"),
      legend.position = "right"
    ) +
    ggplot2::labs(
      x = bquote(Log[2]~Fold~Change~(with~.(conf_pct)*"%"~CI)),
      y = "",
      title = title_txt,
      subtitle = subtitle_txt
    )

  tag_str <- if(!is.null(tag) && nzchar(tag)) paste0("_", tag) else ""
  file_name <- glue::glue("{analyte}_{quant_source}_candidates_{str_quo_name}{tag_str}")

  return(plot_save(
    plot,
    data = if ("tidyproteomics" %in% class(data)) data else NULL,
    file_name = file_name,
    destination = destination,
    height = height,
    width = width,
    path = path
  ))
}

#' @rdname plot_candidates
#' @export
plot_forest <- function(...) {
  plot_candidates(...)
}
