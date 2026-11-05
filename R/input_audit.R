#' Audit raw data headers against tidyproteomics platform configuration
#'
#' @description
#' Inspects the column headers of a raw proteomics data file, matches them
#' against the regex patterns defined in the platform configuration, and
#' produces a comprehensive audit table with rules-based recommendations for
#' unmapped or alternative columns (e.g. gene names, descriptions, MaxLFQ intensity).
#'
#' @param file_name Path to the raw data file (txt, tsv, csv, etc.)
#' @param platform The source of the data (e.g., 'MaxQuant', 'ProteomeDiscoverer', 'FragPipe', 'DIA-NN')
#' @param analyte The omics analyte ('proteins', 'peptides')
#' @param config Optional path to a specific config TSV file, or a data.frame of config rows
#'
#' @return A tibble with raw headers, matching status, config mapping, and suggestions
#' @export
#'
#' @importFrom magrittr %>%
#'
#' @examples
#' library(tidyproteomics)
#' \dontrun{
#' audit <- input_audit("proteinGroups.txt", platform = "MaxQuant", analyte = "proteins")
#' print(audit)
#' }
#'
input_audit <- function(
    file_name = NULL,
    platform = NULL,
    analyte = c('proteins', 'peptides'),
    config = NULL
){
  if (is.null(file_name) || !nzchar(file_name)) {
    cli::cli_abort(c("x" = "No file_name provided for input_audit()."))
  }
  if (!file.exists(file_name)) {
    cli::cli_abort(c("x" = glue::glue("File '{file_name}' not found.")))
  }

  analyte <- rlang::arg_match(analyte)

  # 1. Load config table
  tbl_config <- NULL
  if (is.null(config)) {
    local_cfg <- glue::glue("{platform}_{analyte}.tsv")
    if (file.exists(local_cfg)) {
      tbl_config <- utils::read.delim(local_cfg, sep = "	", stringsAsFactors = FALSE)
      cli::cli_alert_info("Using local configuration file: {.file {local_cfg}}")
    } else {
      if (is.null(platform) || !nzchar(platform)) {
        cli::cli_abort(c("x" = "No platform specified for input_audit()."))
      }
      path_to_config <- glue::glue("{platform}_{analyte}") %>% path_to_package_data()
      if (!file.exists(path_to_config)) {
        cli::cli_abort(c("x" = glue::glue("Config file not found for {platform}_{analyte}.")))
      }
      tbl_config <- utils::read.delim(path_to_config, sep = "	", stringsAsFactors = FALSE)
      cli::cli_alert_info("Using package configuration: {.emph {platform}_{analyte}}")
    }
  } else if (is.character(config) && file.exists(config)) {
    tbl_config <- utils::read.delim(config, sep = "	", stringsAsFactors = FALSE)
    cli::cli_alert_info("Using specified configuration file: {.file {config}}")
  } else if (is.data.frame(config)) {
    tbl_config <- config
    cli::cli_alert_info("Using in-memory configuration data.frame")
  } else {
    cli::cli_abort(c("x" = glue::glue("Invalid config provided: '{config}'.")))
  }

  # Ensure expected columns in config
  req_cols <- c("category", "column_defined", "column_import")
  if (!all(req_cols %in% names(tbl_config))) {
    cli::cli_abort(c("x" = "Config table must contain: category, column_defined, column_import."))
  }

  # 2. Read first lines of raw file to extract headers
  first_lines <- readLines(file_name, n = 2, warn = FALSE)
  if (length(first_lines) == 0) {
    cli::cli_abort(c("x" = glue::glue("Raw data file '{file_name}' is empty.")))
  }
  delim <- if (grepl("	", first_lines[1])) "	" else if (grepl(",", first_lines[1])) "," else ";"
  raw_headers <- unlist(strsplit(first_lines[1], delim, fixed = TRUE))
  raw_headers <- trimws(gsub('^["\']|["\']$', "", raw_headers))

  row1_vals <- if (length(first_lines) > 1) {
    vals <- unlist(strsplit(first_lines[2], delim, fixed = TRUE))
    trimws(gsub('^["\']|["\']$', "", vals))
  } else {
    rep("", length(raw_headers))
  }

  # 3. Match headers against config regex
  has_lfq <- any(grepl("^LFQ\\s+intensity\\s", raw_headers, ignore.case = TRUE))
  has_raw_intensity <- any(grepl("^Intensity\\s", raw_headers, ignore.case = TRUE))

  matched_cfg <- tbl_config[!is.na(tbl_config$column_import) & nzchar(tbl_config$column_import), ]

  results <- list()
  for (i in seq_along(raw_headers)) {
    h <- raw_headers[i]
    val <- if (i <= length(row1_vals)) row1_vals[i] else ""

    is_match <- FALSE
    m_cat <- NA_character_
    m_def <- NA_character_
    m_reg <- NA_character_
    sugg <- NA_character_

    for (r in seq_len(nrow(matched_cfg))) {
      pat <- matched_cfg$column_import[r]
      if (grepl(pat, h, perl = TRUE)) {
        is_match <- TRUE
        m_cat <- matched_cfg$category[r]
        m_def <- matched_cfg$column_defined[r]
        m_reg <- pat
        break
      }
    }

    if (is_match) {
      if (m_def == "abundance_raw" && has_lfq && grepl("^Intensity\\s", m_reg)) {
        sugg <- "Mapped to raw Intensity. MaxLFQ columns detected; switch to '^LFQ intensity\\s' if normalized intensity is desired."
      } else {
        sugg <- "Matched"
      }
    } else {
      # Rules-based suggestions
      if (grepl("(?i)^(gene(\\s*names?|\\s*symbols?)?|gn)$", h)) {
        sugg <- "Suggested annotation::gene_name"
      } else if (grepl("(?i)^(protein\\s*names?|descriptions?|protein\\s*description)$", h)) {
        sugg <- "Suggested annotation::description"
      } else if (grepl("(?i)^(fasta\\s*headers?)$", h)) {
        sugg <- "Suggested annotation::fasta_header"
      } else if (grepl("(?i)^lfq\\s*intensity\\s", h)) {
        sugg <- "Candidate quantitative::abundance_raw (MaxLFQ normalized intensity)"
      } else if (grepl("(?i)(majority\\s*protein\\s*ids?|leading\\s*protein|protein\\s*accessions?|accessions?)", h)) {
        sugg <- "Candidate identifier::protein (alternative accession)"
      } else if (grepl("(?i)(potential\\s*contaminant|reverse|only\\s*identified\\s*by\\s*site)", h)) {
        sugg <- "Candidate filter flag"
      } else if (grepl("(?i)(ms/ms\\s*count|psms?|spectral\\s*counts?|unique\\s*peptides?|razor)", h)) {
        sugg <- "Candidate accounting column"
      } else if (grepl("(?i)(entrez|ensembl|uniprot|organism|species|taxonomy)", h)) {
        sugg <- "Candidate annotation column"
      } else {
        sugg <- "Unassigned"
      }
    }

    results[[length(results) + 1]] <- data.frame(
      raw_header = h,
      matched = is_match,
      category = m_cat,
      column_defined = m_def,
      config_regex = m_reg,
      suggestion = sugg,
      stringsAsFactors = FALSE
    )
  }

  audit_tbl <- tibble::as_tibble(do.call(rbind, results))

  # Summary reporting
  total_col <- nrow(audit_tbl)
  matched_count <- sum(audit_tbl$matched)
  unmatched_count <- total_col - matched_count

  cli::cli_div(theme = list(span.emph = list(color = "blue"), span.info = list(color = "green")))
  cli::cli_h1("Input Header Audit Summary")
  cli::cli_alert_info("File: {.file {basename(file_name)}} ({total_col} columns)")
  cli::cli_alert_success("Matched Columns: {.info {matched_count}} / {total_col}")

  # Check key mandatory fields
  id_matched <- audit_tbl %>% dplyr::filter(category == "identifier")
  quant_matched <- audit_tbl %>% dplyr::filter(category == "quantitative")
  annot_matched <- audit_tbl %>% dplyr::filter(category == "annotation")

  if (nrow(id_matched) > 0) {
    cli::cli_alert_success("Identifier: {.emph {id_matched$raw_header[1]}} -> {.info {id_matched$column_defined[1]}}")
  } else {
    cli::cli_alert_danger("Missing required identifier column mapping!")
  }

  if (nrow(quant_matched) > 0) {
    cli::cli_alert_success("Quantitative: {.emph {length(unique(quant_matched$raw_header))}} columns matched to {.info {unique(quant_matched$column_defined)}}")
  } else {
    cli::cli_alert_danger("Missing required quantitative column mapping!")
  }

  if (nrow(annot_matched) > 0) {
    cli::cli_alert_success("Annotations: {.emph {paste(unique(annot_matched$column_defined), collapse = ', ')}}")
  }

  # Actionable suggestions
  actionable <- audit_tbl %>% dplyr::filter(!matched & suggestion != "Unassigned")
  if (nrow(actionable) > 0) {
    cli::cli_h2("Actionable Suggestions")
    for (s in unique(actionable$suggestion)) {
      cols_s <- actionable$raw_header[actionable$suggestion == s]
      cli::cli_alert_info("{s}: {.emph {paste(head(cols_s, 5), collapse = ', ')}}")
    }
  }

  return(audit_tbl)
}
