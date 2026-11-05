#' Plot Biological Signatures or Program Dichotomies
#'
#' @description
#' `plot_signatures()` visualizes directional expression shifts across distinct
#' biological modules, pathways, or opposing signatures (e.g., comparing upregulated
#' quality control machinery against downregulated cell-cycle/growth machinery)
#' as a faceted, diverging lollipop or bar chart.
#'
#' @param data a tidyproteomics data object (or an exported expression table)
#' @param ... two sample comparison expression (e.g. `knockdown/control`)
#' @param signatures a named list of character vectors of identifiers/gene names
#'   representing distinct biological signatures or programs (e.g. `list("UPR" = c(...), "CellCycle" = c(...))`).
#' @param programs alias for `signatures` for backward compatibility.
#' @param set_a optional convenience parameter for the first program gene set (used if `signatures` is NULL).
#' @param set_b optional convenience parameter for the second program gene set (used if `signatures` is NULL).
#' @param label_a character label for `set_a`. Default is "Signature A".
#' @param label_b character label for `set_b`. Default is "Signature B".
#' @param labels_column a character defining the column name for feature labeling. Default is "gene_name".
#' @param log2fc_column a character defining the column name of log2 foldchange values. Default is "log2_foldchange".
#' @param significance_column a character defining the column name of significance values. Default is "adj_p_value".
#' @param significance_max a numeric defining the significance threshold. Default is 0.05.
#' @param plot_type character defining geometry type, either "lollipop" (default) or "bar".
#' @param show_title boolean or character string for plot title. Default is TRUE.
#' @param tag optional character string to append to the plot file name.
#' @param destination character string indicating output type: "plot", "save", "png", "svg", "tiff", "jpeg". Default is "plot".
#' @param height numeric plot height in inches. Default is 6.
#' @param width numeric plot width in inches. Default is 8.
#' @param path optional path to save destination plot.
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
#'   plot_signatures(
#'     knockdown/control,
#'     signatures = list(
#'       "Ribosomal" = c("RPL26", "RPS27A"),
#'       "Cell Cycle" = c("CKS2", "RACGAP1")
#'     ),
#'     labels_column = "gene_name"
#'   )
#'
plot_signatures <- function(
    data = NULL,
    ...,
    signatures = NULL,
    programs = NULL,
    set_a = NULL,
    set_b = NULL,
    label_a = "Signature A",
    label_b = "Signature B",
    labels_column = "gene_name",
    log2fc_column = "log2_foldchange",
    significance_column = "adj_p_value",
    significance_max = 0.05,
    plot_type = c("lollipop", "bar"),
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
  signature <- NULL
  is_sig <- NULL

  plot_type <- match.arg(plot_type)

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

  # resolve signatures alias with programs
  if (is.null(signatures) && !is.null(programs)) {
    signatures <- programs
  }

  # resolve signatures specification
  if (is.null(signatures)) {
    if (!is.null(set_a) && !is.null(set_b)) {
      signatures <- list(set_a, set_b)
      names(signatures) <- c(label_a, label_b)
    } else {
      cli::cli_abort("Either `signatures` (as a named list) or both `set_a` and `set_b` must be supplied.")
    }
  }

  if (!is.list(signatures) || is.null(names(signatures)) || any(!nzchar(names(signatures)))) {
    cli::cli_abort("`signatures` must be a named list of character vectors (e.g. list('Signature A' = c(...), 'Signature B' = c(...)))")
  }

  # filter and assemble signature subsets
  sig_dfs <- lapply(names(signatures), function(s_name) {
    targets <- as.character(signatures[[s_name]])
    sub <- table %>%
      dplyr::filter(
        (!is.na(.data[[label_var]]) & .data[[label_var]] %in% targets) |
          (!is.na(.data[[id_col]]) & .data[[id_col]] %in% targets)
      )
    if (nrow(sub) > 0) {
      sub$signature <- s_name
    }
    sub
  })

  sig_table <- dplyr::bind_rows(sig_dfs)

  if (nrow(sig_table) == 0) {
    cli::cli_alert_warning("No matching proteins found in data for the supplied signatures.")
    return(if ("tidyproteomics" %in% class(data)) data else NULL)
  }

  # clean labels and order factors
  raw_labels <- sig_table[[label_var]]
  labels_clean <- ifelse(is.na(raw_labels) | trimws(as.character(raw_labels)) == "",
                         sig_table[[id_col]],
                         as.character(raw_labels))
  labels_clean <- stringr::str_extract(labels_clean, "^[^;]+(;\\s*[^;]+)?")
  sig_table$plot_label <- labels_clean
  sig_table$is_sig <- ifelse(sig_table[[significance_column]] <= significance_max, "TRUE", "FALSE")

  # order signatures as supplied and order proteins within signatures by effect size
  sig_table$signature <- factor(sig_table$signature, levels = names(signatures))
  sig_table <- sig_table %>%
    dplyr::arrange(signature, .data[[log2fc_column]])

  # factor plot_label preserving ranked order
  unique_labels <- unique(sig_table$plot_label)
  sig_table$plot_label <- factor(sig_table$plot_label, levels = unique_labels)

  # build palette
  num_sigs <- length(names(signatures))
  palette_colors <- theme_palette(max(num_sigs, 2))[1:num_sigs]
  names(palette_colors) <- names(signatures)

  # titles
  title_txt <- glue::glue("{stringr::str_to_title(gsub('_', ' ', str_quo_name))} — Biological Signatures")
  if (is.character(show_title)) {
    title_txt <- show_title
  } else if (isFALSE(show_title)) {
    title_txt <- NULL
  }

  subtitle_txt <- glue::glue("{nrow(sig_table)} proteins across {num_sigs} signatures | {significance_column} <= {significance_max}")
  if (isFALSE(show_title)) {
    subtitle_txt <- NULL
  }

  # ggplot construction
  plot <- ggplot2::ggplot(
    sig_table,
    ggplot2::aes(
      x = .data[[log2fc_column]],
      y = plot_label,
      color = signature,
      fill = signature
    )
  ) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "grey50", alpha = 0.6)

  if (plot_type == "lollipop") {
    plot <- plot +
      ggplot2::geom_segment(
        ggplot2::aes(x = 0, xend = .data[[log2fc_column]], y = plot_label, yend = plot_label),
        linewidth = 0.8,
        alpha = 0.7
      ) +
      ggplot2::geom_point(
        ggplot2::aes(shape = is_sig),
        size = 3
      ) +
      ggplot2::scale_shape_manual(
        values = c("TRUE" = 16, "FALSE" = 1),
        name = glue::glue("Sig ({significance_max})"),
        drop = FALSE
      )
  } else {
    plot <- plot +
      ggplot2::geom_col(
        width = 0.6,
        alpha = 0.85
      )
  }

  plot <- plot +
    ggplot2::facet_grid(
      signature ~ .,
      scales = "free_y",
      space = "free_y",
      labeller = ggplot2::labeller(signature = ggplot2::label_wrap_gen(width = 22))
    ) +
    ggplot2::scale_color_manual(values = palette_colors) +
    ggplot2::scale_fill_manual(values = palette_colors) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      strip.text = ggplot2::element_text(face = "bold", size = 8.5),
      strip.background = ggplot2::element_rect(fill = "grey95"),
      legend.position = "right"
    ) +
    ggplot2::labs(
      x = expression(Log[2]~Fold~Change),
      y = "",
      title = title_txt,
      subtitle = subtitle_txt
    )

  tag_str <- if(!is.null(tag) && nzchar(tag)) paste0("_", tag) else ""
  file_name <- glue::glue("{analyte}_{quant_source}_signatures_{str_quo_name}{tag_str}")

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

#' @rdname plot_signatures
#' @export
plot_programs <- function(...) {
  plot_signatures(...)
}
