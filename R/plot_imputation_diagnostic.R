#' Plot Imputation Diagnostic Across Conditions
#'
#' @description
#' `plot_imputation_diagnostic()` provides a targeted quality diagnostic visual for
#' missing-value imputation. It visualizes individual replicate abundances across two
#' conditions for candidate proteins, distinguishing observed values from those imputed
#' (e.g. at the limit of detection). By default, it prioritizes top candidates that possess
#' imputed values, allowing rapid assessment of whether fold changes are continuous or
#' detection-limited on/off artifacts.
#'
#' @param data a tidyproteomics data object
#' @param ... two sample comparison expression (e.g. `knockdown/control`)
#' @param candidates a character vector of protein identifiers or gene names to plot. If NULL,
#'   defaults to `proteins` or automatically selects the top candidates with imputed values.
#' @param proteins alias for `candidates` for backward compatibility.
#' @param labels_column a character defining the column name for facet titles (default is `"gene_name"`).
#' @param top_n integer defining the number of proteins to plot when `candidates` is NULL. Default is 10.
#' @param log2fc_column a character defining the column name of log2 foldchange values. Default is `"log2_foldchange"`.
#' @param significance_column a character defining the column name of significance values. Default is `"adj_p_value"`.
#' @param show_title boolean or character string for plot title. Default is TRUE.
#' @param tag optional character string to append to the plot file name.
#' @param destination character string indicating output type: `"plot"`, `"save"`, `"png"`, `"svg"`, `"tiff"`, `"jpeg"`. Default is `"plot"`.
#' @param height numeric plot height in inches. Default is 6.
#' @param width numeric plot width in inches. Default is 8.
#' @param path optional path to save destination plot.
#' @param ... passthrough arguments.
#'
#' @return a ggplot2 object or tidyproteomics data object depending on `destination`.
#' @export
#'
#' @examples
#' library(dplyr, warn.conflicts = FALSE)
#' library(tidyproteomics)
#'
#' # subset for fast execution
#' hela_proteins <- hela_proteins %>% subset(description %like% "Ribosome")
#'
#' hela_proteins %>%
#'   expression(knockdown/control) %>%
#'   plot_imputation_diagnostic(knockdown/control, top_n = 6)
#'
plot_imputation_diagnostic <- function(
    data = NULL,
    ...,
    candidates = NULL,
    proteins = NULL,
    labels_column = "gene_name",
    top_n = 10,
    log2fc_min = 1,
    log2fc_column = "log2_foldchange",
    significance_max = 0.05,
    significance_column = "adj_p_value",
    regulation = c("both", "up", "down"),
    show_title = TRUE,
    tag = NULL,
    destination = "plot",
    height = 6,
    width = 8,
    path = NULL
){
  check_data(data)

  str_quo <- tidyproteomics_quo(...)
  if(is.null(str_quo)) {
    cli::cli_alert_warning("No comparison provided")
    cli::cli_alert_info("Try one of ({paste(names(data$analysis), collapse = ', ')})")
    return(data)
  }

  str_quo_name <- tidyproteomics_quo_name(...)
  tag_str <- if(!is.null(tag) && nzchar(tag)) paste0("_", tag) else ""
  file_name <- glue::glue("{data$analyte}_{data$quantitative_source}_imputation_diagnostic_{str_quo_name}{tag_str}")

  title_arg <- show_title
  if (isTRUE(show_title)) {
    base_title <- str_quo_name %>%
      stringr::str_replace("-", " / ") %>%
      stringr::str_replace_all("_", " ") %>%
      stringr::str_to_title()
    title_arg <- as.character(glue::glue("{base_title} — Imputation Diagnostic"))
  }

  plot_expression(
    data = data,
    ...,
    tag = tag,
    log2fc_min = log2fc_min,
    log2fc_column = log2fc_column,
    significance_max = significance_max,
    significance_column = significance_column,
    top_n = top_n,
    regulation = regulation,
    labels_column = labels_column,
    candidates = candidates,
    proteins = proteins,
    imputed_only = TRUE,
    show_title = title_arg,
    destination = destination,
    height = height,
    width = width,
    path = path,
    file_name = file_name
  )
}
