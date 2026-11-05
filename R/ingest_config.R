#' Edit or customize a platform ingestion configuration table
#'
#' @description
#' Facilitates viewing, adding, adjusting, or auto-configuring regex patterns and
#' annotation mappings in a platform configuration TSV file. Works seamlessly with
#' `export_config()` and `input_audit()`.
#'
#' @param platform Platform source (e.g. 'MaxQuant', 'ProteomeDiscoverer', 'FragPipe', 'DIA-NN')
#' @param analyte Omics analyte ('proteins', 'peptides')
#' @param config_file Optional path to a specific config TSV file. If NULL, targets `<platform>_<analyte>.tsv` in current directory.
#' @param action Action to execute: 'view', 'set', 'add', 'remove', 'auto_fix'
#' @param category Category for the entry ('identifier', 'quantitative', 'accounting', 'annotation', 'sample', 'impute')
#' @param column_defined Variable name in tidyproteomics (e.g. 'gene_name', 'description', 'abundance_raw', 'protein')
#' @param column_import Regex pattern to match the raw file column header (e.g. '^Gene names$', '^LFQ intensity\\s')
#' @param pattern_extract Optional regex pattern to extract substring
#' @param pattern_remove Optional regex pattern to remove
#' @param pattern_split Optional regex pattern to split (e.g. '\\;')
#' @param pivot Optional boolean (TRUE for quantitative/accounting sample-pivoted columns)
#' @param file_name Optional raw data file path required when action = 'auto_fix'
#' @param use_lfq Optional logical when action = 'auto_fix', whether to prioritize MaxLFQ intensity (default: TRUE if LFQ present)
#' @param save Logical, whether to write the updated configuration to disk (default: TRUE)
#' @param output_file Output TSV path. Defaults to `config_file` or `<platform>_<analyte>.tsv`.
#'
#' @return A tibble of the modified configuration table (invisibly)
#' @export
#'
#' @importFrom magrittr %>%
#'
#' @examples
#' library(tidyproteomics)
#' \dontrun{
#' # View current config
#' ingest_config("MaxQuant", "proteins", action = "view")
#'
#' # Add an annotation
#' ingest_config("MaxQuant", "proteins", action = "add",
#'               category = "annotation", column_defined = "gene_name",
#'               column_import = "^Gene names$")
#'
#' # Update quantitative layer to MaxLFQ
#' ingest_config("MaxQuant", "proteins", action = "set",
#'               column_defined = "abundance_raw", column_import = "^LFQ intensity\\s")
#' }
#'
ingest_config <- function(
    platform = NULL,
    analyte = c('proteins', 'peptides'),
    config_file = NULL,
    action = c('view', 'set', 'add', 'remove', 'auto_fix'),
    category = NULL,
    column_defined = NULL,
    column_import = NULL,
    pattern_extract = NULL,
    pattern_remove = NULL,
    pattern_split = NULL,
    pivot = NULL,
    file_name = NULL,
    use_lfq = NULL,
    save = TRUE,
    output_file = NULL
){
  action <- rlang::arg_match(action)
  analyte <- rlang::arg_match(analyte)

  # 1. Establish the configuration file
  target_file <- if (!is.null(config_file) && nzchar(config_file)) {
    config_file
  } else {
    glue::glue("{platform}_{analyte}.tsv")
  }

  tbl_config <- NULL
  if (file.exists(target_file)) {
    tbl_config <- utils::read.delim(target_file, stringsAsFactors = FALSE)
  } else {
    if (is.null(platform) || !nzchar(platform)) {
      cli::cli_abort(c("x" = "No platform specified and '{target_file}' does not exist."))
    }
    path_to_config <- glue::glue("{platform}_{analyte}") %>% path_to_package_data()
    if (!file.exists(path_to_config)) {
      cli::cli_abort(c("x" = glue::glue("Config file not found for {platform}_{analyte}.")))
    }
    file.copy(path_to_config, target_file)
    cli::cli_alert_success("Exported base configuration to {.file {target_file}}")
    tbl_config <- utils::read.delim(target_file, stringsAsFactors = FALSE)
  }

  # Ensure expected schema
  req_cols <- c("category", "column_defined", "column_import", "pattern_extract", "pattern_remove", "pattern_split", "pivot")
  for (col in req_cols) {
    if (!col %in% names(tbl_config)) {
      tbl_config[[col]] <- NA
    }
  }

  # 2. Execute requested action
  if (action == "view") {
    cli::cli_h1("Configuration: {basename(target_file)}")
    print(tibble::as_tibble(tbl_config))
    return(invisible(tibble::as_tibble(tbl_config)))
  }

  if (action == "add") {
    if (is.null(column_defined) || !nzchar(column_defined)) {
      cli::cli_abort(c("x" = "column_defined must be specified to add a configuration row."))
    }
    if (is.null(column_import) || !nzchar(column_import)) {
      cli::cli_abort(c("x" = "column_import regex must be specified to add a configuration row."))
    }
    cat_val <- if (!is.null(category) && nzchar(category)) category else "annotation"

    existing_idx <- which(tbl_config$column_defined == column_defined)
    if (length(existing_idx) > 0) {
      cli::cli_alert_info("Column {.emph {column_defined}} already exists in config. Updating existing row.")
      tbl_config$category[existing_idx] <- cat_val
      tbl_config$column_import[existing_idx] <- column_import
      if (!is.null(pattern_extract)) tbl_config$pattern_extract[existing_idx] <- pattern_extract
      if (!is.null(pattern_remove)) tbl_config$pattern_remove[existing_idx] <- pattern_remove
      if (!is.null(pattern_split)) tbl_config$pattern_split[existing_idx] <- pattern_split
      if (!is.null(pivot)) tbl_config$pivot[existing_idx] <- pivot
    } else {
      new_row <- data.frame(
        category = cat_val,
        column_defined = column_defined,
        column_import = column_import,
        pattern_extract = if (!is.null(pattern_extract)) pattern_extract else NA_character_,
        pattern_remove = if (!is.null(pattern_remove)) pattern_remove else NA_character_,
        pattern_split = if (!is.null(pattern_split)) pattern_split else NA_character_,
        pivot = if (!is.null(pivot)) pivot else NA,
        stringsAsFactors = FALSE
      )
      tbl_config <- rbind(tbl_config, new_row)
      cli::cli_alert_success("Added {.info {cat_val}}::{.emph {column_defined}} with pattern {.code {column_import}}")
    }
  }

  if (action == "set") {
    if (is.null(column_defined) || !nzchar(column_defined)) {
      cli::cli_abort(c("x" = "column_defined must be specified to update a configuration row."))
    }
    idx <- which(tbl_config$column_defined == column_defined)
    if (length(idx) == 0) {
      cli::cli_abort(c("x" = glue::glue("Column '{column_defined}' not found in config. Use action = 'add' to create it.")))
    }
    if (!is.null(category)) tbl_config$category[idx] <- category
    if (!is.null(column_import)) tbl_config$column_import[idx] <- column_import
    if (!is.null(pattern_extract)) tbl_config$pattern_extract[idx] <- pattern_extract
    if (!is.null(pattern_remove)) tbl_config$pattern_remove[idx] <- pattern_remove
    if (!is.null(pattern_split)) tbl_config$pattern_split[idx] <- pattern_split
    if (!is.null(pivot)) tbl_config$pivot[idx] <- pivot

    # Smart helper: If setting abundance_raw to LFQ intensity, ensure sample and sample_file regexes match LFQ
    if (column_defined == "abundance_raw" && !is.null(column_import) && grepl("LFQ", column_import, ignore.case = TRUE)) {
      sample_idx <- which(tbl_config$category == "sample" & tbl_config$column_defined %in% c("sample", "sample_file"))
      if (length(sample_idx) > 0) {
        updated_extract <- "(?<=(LFQ\\sintensity|Intensity|Peptides|Unique\\speptides)\\s).+"
        for (s_idx in sample_idx) {
          curr_extract <- tbl_config$pattern_extract[s_idx]
          if (!is.na(curr_extract) && !grepl("LFQ", curr_extract, ignore.case = TRUE)) {
            tbl_config$pattern_extract[s_idx] <- updated_extract
          }
        }
        cli::cli_alert_info("Synchronized sample/sample_file pattern_extract to: {.code {updated_extract}}")
      }
    }
    cli::cli_alert_success("Updated {.emph {column_defined}} in configuration.")
  }

  if (action == "remove") {
    if (is.null(column_defined) || !nzchar(column_defined)) {
      cli::cli_abort(c("x" = "column_defined must be specified to remove a configuration row."))
    }
    idx <- which(tbl_config$column_defined == column_defined)
    if (length(idx) == 0) {
      cli::cli_alert_warning("Column '{column_defined}' was not found in config.")
    } else {
      tbl_config <- tbl_config[-idx, , drop = FALSE]
      cli::cli_alert_success("Removed {.emph {column_defined}} from configuration.")
    }
  }

  if (action == "auto_fix") {
    if (is.null(file_name) || !file.exists(file_name)) {
      cli::cli_abort(c("x" = "action = 'auto_fix' requires a valid file_name to inspect raw headers."))
    }
    audit_res <- input_audit(file_name = file_name, platform = platform, analyte = analyte, config = tbl_config)
    raw_headers <- audit_res$raw_header

    # Auto-add gene_name if present and missing
    if (!"gene_name" %in% tbl_config$column_defined) {
      g_col <- grep("(?i)^(gene(\\s*names?|\\s*symbols?)?|gn)$", raw_headers, value = TRUE)
      if (length(g_col) > 0) {
        tbl_config <- rbind(tbl_config, data.frame(
          category = "annotation",
          column_defined = "gene_name",
          column_import = glue::glue("^{g_col[1]}$"),
          pattern_extract = NA_character_, pattern_remove = NA_character_, pattern_split = NA_character_, pivot = NA,
          stringsAsFactors = FALSE
        ))
        cli::cli_alert_success("Auto-configured annotation: {.emph gene_name} -> {.code ^{g_col[1]}$}")
      }
    }

    # Auto-add description if present and missing
    if (!"description" %in% tbl_config$column_defined) {
      d_col <- grep("(?i)^(protein\\s*names?|descriptions?|protein\\s*description)$", raw_headers, value = TRUE)
      if (length(d_col) > 0) {
        tbl_config <- rbind(tbl_config, data.frame(
          category = "annotation",
          column_defined = "description",
          column_import = glue::glue("^{d_col[1]}$"),
          pattern_extract = NA_character_, pattern_remove = NA_character_, pattern_split = NA_character_, pivot = NA,
          stringsAsFactors = FALSE
        ))
        cli::cli_alert_success("Auto-configured annotation: {.emph description} -> {.code ^{d_col[1]}$}")
      }
    }

    # Auto-add fasta_header if present and missing
    if (!"fasta_header" %in% tbl_config$column_defined) {
      f_col <- grep("(?i)^(fasta\\s*headers?)$", raw_headers, value = TRUE)
      if (length(f_col) > 0) {
        tbl_config <- rbind(tbl_config, data.frame(
          category = "annotation",
          column_defined = "fasta_header",
          column_import = glue::glue("^{f_col[1]}$"),
          pattern_extract = NA_character_, pattern_remove = NA_character_, pattern_split = NA_character_, pivot = NA,
          stringsAsFactors = FALSE
        ))
        cli::cli_alert_success("Auto-configured annotation: {.emph fasta_header} -> {.code ^{f_col[1]}$}")
      }
    }

    # Auto-switch to MaxLFQ if requested or detected
    want_lfq <- if (!is.null(use_lfq)) isTRUE(use_lfq) else any(grepl("^LFQ\\s+intensity\\s", raw_headers, ignore.case = TRUE))
    if (want_lfq && any(grepl("^LFQ\\s+intensity\\s", raw_headers, ignore.case = TRUE))) {
      q_idx <- which(tbl_config$column_defined == "abundance_raw")
      if (length(q_idx) > 0 && !grepl("LFQ", tbl_config$column_import[q_idx], ignore.case = TRUE)) {
        tbl_config$column_import[q_idx] <- "^LFQ intensity\\s"
        cli::cli_alert_success("Switched quantitative metric: {.emph abundance_raw} -> {.code ^LFQ intensity\\s}")
        sample_idx <- which(tbl_config$category == "sample" & tbl_config$column_defined %in% c("sample", "sample_file"))
        if (length(sample_idx) > 0) {
          updated_extract <- "(?<=(LFQ\\sintensity|Intensity|Peptides|Unique\\speptides)\\s).+"
          for (s_idx in sample_idx) {
            if (!is.na(tbl_config$pattern_extract[s_idx]) && nzchar(tbl_config$pattern_extract[s_idx])) {
              tbl_config$pattern_extract[s_idx] <- updated_extract
            }
          }
          cli::cli_alert_info("Synchronized sample/sample_file pattern_extract to: {.code {updated_extract}}")
        }
      }
    }
  }

  # 3. Save to disk if requested
  out_path <- if (!is.null(output_file) && nzchar(output_file)) {
    output_file
  } else {
    target_file
  }

  if (isTRUE(save)) {
    utils::write.table(tbl_config, file = out_path, sep = "	", row.names = FALSE, quote = FALSE, na = "")
    cli::cli_alert_success("Configuration written to {.file {out_path}}")
  }

  return(invisible(tibble::as_tibble(tbl_config)))
}
