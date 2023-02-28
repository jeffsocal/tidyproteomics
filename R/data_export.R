#' Export the quantitative data from an tidyproteomics data-object
#'
#' @description
#' `export_quant()` returns the main quantitative data object as a tibble with
#' _identifier_ as the designation for the measured observation.
#'
#' @param data tidyproteomics data object
#' @param file_name character string vector
#' @param raw_data a boolean
#' @param normalized a boolean
#' @param scaled a boolean
#'
#' @return a tibble
#' @export
#'
#' @examples
#' library(dplyr, warn.conflicts = FALSE)
#' library(tidyproteomics)
#' hela_proteins %>%
#'    normalize(.method = "loess") %>%
#'    export_quant(file_name = "hela_quant_data.xlsx", normalized = "loess")
#'
export_quant <- function(
    data       = NULL,
    file_name  = NULL,
    raw_data   = TRUE,
    normalized = FALSE,
    scaled     = c('none', 'between', 'proportion')
){

  # visible bindings
  abundance <- NULL
  abundance_raw <- NULL
  abundance_norm <- NULL
  identifier <- NULL
  abundance_scaled <- NULL

  check_data(data)
  scaled <- rlang::arg_match(scaled)
  analyte <- data$analyte
  analyte <- rlang::arg_match(analyte, c('proteins','peptides'))

  cli::cli_div(theme = list(span.emph = list(color = "#ff4500")))

  data_quant <- extract(data, 'raw') %>%
    dplyr::rename(abundance_raw = abundance) %>%
    dplyr::filter(!is.na(abundance_raw))

  if(normalized != FALSE) {

    data_norm <- extract(data, normalized)

    data_quant <- data_quant %>%
      dplyr::full_join(data_norm,
                       by=c('identifier', 'sample', 'replicate')) %>%
      dplyr::rename(abundance_norm = abundance) %>%
      dplyr::filter(!is.na(abundance_norm))

  }
  #
  # pull in the normalized values
  #
  if(scaled != FALSE){
    if(normalized != FALSE){
      data_quant <- data_quant %>%
        dplyr::mutate(abundance_scaled = abundance_norm)
    } else {
      data_quant <- data_quant %>%
        dplyr::mutate(abundance_scaled = abundance_raw)
    }

    if(scaled == 'between'){
      data_quant <- data_quant %>%
        dplyr::group_by(identifier) %>%
        dplyr::mutate(abundance_scaled = abundance_scaled / sum(abundance_scaled) * 100) %>%
        dplyr::ungroup()
    } else {
      data_quant <- data_quant %>%
        dplyr::group_by(sample, replicate) %>%
        dplyr::mutate(abundance_scaled = abundance_scaled / sum(abundance_scaled) * 100) %>%
        dplyr::ungroup()
    }

    data_quant <- data_quant %>% dplyr::rename_at(dplyr::vars(tidyselect::starts_with("abundance_scaled")), ~paste('abundance_scaled', scaled, sep="-"))
  }

  if(normalized != FALSE){
    data_quant <- data_quant %>% dplyr::rename_at(dplyr::vars(tidyselect::starts_with("abundance_norm")), ~paste('abundance', normalized, sep="_"))
  }

  col_abn <- names(data_quant)[which(grepl("abundance", names(data_quant)))]

  tbl_out <- data_quant %>%
    tidyr::unite(sample, sample, replicate) %>%
    tidyr::pivot_wider(
      names_from = 'sample',
      values_from = tidyselect::all_of(col_abn)
    ) %>%
    munge_identifier("separate", data$identifier) %>%
    dplyr::left_join(
      data$annotations %>%
        tidyr::pivot_wider(names_from = 'term', values_from = 'annotation'),
      by = data$identifier
    )

  if(is.null(file_name)) {return(tbl_out)}
  tbl_out %>% write_local(file_name)
}

#' Export the quantitative data from an tidyproteomics data-object
#'
#' @description
#' `export()` returns the main quantitative data object as a tibble with
#' _identifier_ as the designation for the measured observation.
#'
#' @param data tidyproteomics data object
#' @param ... two sample comparison e.g. experimental/control
#' @param .analysis a character string for the specific analysis to export
#' @param .term a character string of the term from an enrichment analysis
#' @param .append a character string of the term to append to the output
#' @param .file_name a character string for file to write to, format implied from string ('.rds', '.xlsx', '.csv', '.tsv')
#'
#' @return a tibble
#' @export
#'
#' @examples
#' library(dplyr, warn.conflicts = FALSE)
#' library(tidyproteomics)
#' hela_proteins %>%
#'    expression(knockdown/control) %>%
#'    export_analysis(knockdown/control,
#'                    .analysis = "expression")
#'
#' hela_proteins %>%
#'    export_analysis(.analysis = "counts")
#'
export_analysis <- function(
    data = NULL,
    ...,
    .analysis = NULL,
    .term = NULL,
    .append = NULL,
    .file_name = NULL
){

  # visible bindings
  term <- NULL
  annotation <- NULL

  check_data(data)

  if(.analysis == 'counts'){
    tbl_out <- data %>% analysis_counts()
  } else {

    str_quo <- tidyproteomics_quo(...)
    if(is.null(str_quo)) {
      cli::cli_alert_warning("No comparison provided")
      cli::cli_alert_info("Try one of ({paste(names(data$analysis), collapse = ', ')})")
      return(data)
    }
    experiment <- str_quo[['variable']]
    control <- str_quo[['value']]
    if(str_quo['operator'] != "/") {
      cli::cli_div(theme = list(span.emph = list(color = "#ff4500")))
      cli::cli_abort("Comparison operator must be {.emph \"/\"} (e.g. {experiment}{.emph /}{control})")
    }
    check_samples(data, experiment, control)

    set_expression <- glue::glue("{experiment}/{control}")
    if(is.null(data$analysis[[set_expression]])){
      cli::cli_abort(c("No analysis data for `{experiment}` `{control}`",
                       "i" = "  suggestion: run expression({experiment}/{control})"))
    }


    analyses <- c(data$analysis[[set_expression]], 'counts')
    .analysis <- rlang::arg_match(.analysis, names(analyses))
    cli::cli_div(theme = list(span.emph = list(color = "#ff4500")))

    tbl_out <- analyses[[.analysis]]
    if(.analysis == 'enrichment') {
      .term <- rlang::arg_match(.term, names(tbl_out))
      cli::cli_div(theme = list(span.emph = list(color = "#ff4500")))
      tbl_out <- tbl_out[[.term]]
    }

    if(!is.null(data$annotations)){
      if(length(intersect(data$identifier, colnames(tbl_out))) == length(data$identifier)) {
        tbl_out <- tbl_out %>%
          dplyr::left_join(
            data$annotations %>%
              tidyr::pivot_wider(names_from = 'term', values_from = 'annotation'),
            by = data$identifier
          )
      }
    }
  }

  if(!is.null(.append)) {

    .append <- rlang::arg_match(.append, get_annotation_terms(data))

    tbl_out <- tbl_out %>%
      dplyr::left_join(
        data$annotation %>%
          dplyr::filter(term %in% c(.term, .append)) %>%
          tidyr::pivot_wider(names_from = 'term', values_from = 'annotation') %>%
          dplyr::select(dplyr::matches(paste0("^", paste(c(.term, .append), collapse = "$|^"), "$"))) %>%
          tidyr::separate_rows(!!.term, sep = ";") %>%
          dplyr::rename(annotation = !!.term) %>%
          dplyr::rename(append = !!.append) %>%
          dplyr::filter(!is.na(annotation)) %>%
          dplyr::group_by(annotation) %>%
          dplyr::summarize(
            append_group = paste(append, collapse = ", "),
            .groups = 'drop'
          ), by = 'annotation')

    colnames(tbl_out)[which(colnames(tbl_out) == 'append_group')] <- paste0(.append, "s")
  }

  if(is.null(.file_name)) {return(tbl_out)}
  tbl_out %>% write_local(.file_name)
}


#' Store data locally
#'
#' @description
#' `save_local()` will save the tidyproteomics data-object in the local project,
#' based on the given type in the directory ./data/ as either proteins.rds or
#' peptides.rds. This is a _passthrough_ function as it returns the original
#' tidyproteomics data-object.
#'
#' @param data tidyproteomics data object
#'
#' @return tidyproteomics data object
#' @export
#'
save_local <- function(
    data = NULL
) {
  check_data(data)
  if(!dir.exists("./data/rds")) {dir.create("./data/rds", recursive = T)}

  cli::cli_process_start("Saving {data$analyte} to local dir ./data/rds")
  file_name <- paste0(data$analyte, '.rds')
  saveRDS(data, paste0('./data/rds/', file_name))
}

#' Helper functio to write data table locally
#'
#' @description
#' `write_local()` will save the data table in the local project,
#'
#' @param table a tibble
#' @param file_name a tibble
#'
#' @return tidyproteomics data object
#'
write_local <- function(
    table = NULL,
    file_name = NULL
){

  if(is.null(file_name)) {cli::cli_abort("file_name not provided")}

  as <- stringr::str_extract(file_name, "\\..*?$");
  as <- rlang::arg_match(as, c('.rds', '.xlsx', '.csv', '.tsv'))
  if(as == '.xlsx'){
    table %>% writexl::write_xlsx(file_name)
  } else if(as == '.csv') {
    table %>% readr::write_csv(file_name)
  } else if(as == '.tsv') {
    table %>% readr::write_tsv(file_name)
  } else {
    table %>% saveRDS(file_name)
  }
}
