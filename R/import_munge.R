

#' A helper function for importing peptide table data
#'
#' @param tbl_data a table of imported data
#' @param tbl_config a table of config values
#'
#' @return a tibble
#'
import_rename <- function(
    tbl_data = NULL,
    tbl_config = NULL
){
  if(nrow(tbl_config) == 0) return(tbl_data)
  cols_extracted <- names(tbl_data)
  for(i in row.names(tbl_config)){
    pattern <- tbl_config[i,]$pattern_extract
    col_imp <- tbl_config[i,]$column_import
    col_def <- tbl_config[i,]$column_defined
    w <- which(grepl(col_imp, cols_extracted))
    colnames(tbl_data)[w] <- col_def
  }
  return(tbl_data)
}

#' A helper function for importing peptide table data
#'
#' @param tbl_data a table of imported data
#' @param tbl_config a table of config values
#' @param remove as boolean to determine if the extracted column name should change or copy to a new, retaining the old
#'
#' @return a tibble
#'
import_extract <- function(
    tbl_data = NULL,
    tbl_config = NULL,
    remove = FALSE
){

  # visible bindings
  pattern_extract <- NULL
  column_import <- NULL
  category <- NULL

  tbl_config <- tbl_config %>%
    dplyr::filter(!is.na(pattern_extract)) %>%
    dplyr::filter(!is.na(column_import)) %>%
    dplyr::filter(category != 'impute') %>%
    as.data.frame()

  if(nrow(tbl_config) == 0) return(tbl_data)
  for(i in row.names(tbl_config)){
    pattern <- tbl_config[i,]$pattern_extract
    col_imp <- tbl_config[i,]$column_import
    col_def <- tbl_config[i,]$column_defined

    vect_extracted <- tbl_data %>%
      dplyr::select(dplyr::matches(paste0("^", col_imp, "$"))) %>% unlist() %>%
      stringr::str_extract(pattern)

    if(length(vect_extracted) != nrow(tbl_data)) {
      cli::cli_abort("... issue extracting {.emph {pattern}} from {.emph {col_imp} ~ {col_def}}")
    }
    cli::cli_alert_info("... extracting {.emph {pattern}} from {.emph {col_imp} ~ {col_def}}")

    tbl_data[,col_def] <- vect_extracted
    if(remove == TRUE) {
      col <- which(grepl(col_imp, names(tbl_data)))
      tbl_data[,col] <- NULL
    }
  }
  return(tbl_data)
}

#' A helper function for importing peptide table data
#'
#' @param tbl_data a table of imported data
#' @param tbl_config a table of config values
#'
#' @return a tibble
#'
import_split <- function(
    tbl_data = NULL,
    tbl_config = NULL
){

  # visible bindings
  pattern_split <- NULL

  # create the protein_cluster column anyways
  tbl_data <- tbl_data %>% dplyr::mutate(protein_cluster = paste0('pg', dplyr::row_number()))

  tbl_config <- tbl_config %>% dplyr::filter(!is.na(pattern_split)) %>% as.data.frame()
  if(nrow(tbl_config) == 0) { return(tbl_data) }
  for(i in row.names(tbl_config)){

    tbl_cols <- colnames(tbl_data)

    pattern <- tbl_config[i,]$pattern_split
    col_imp <- tbl_config[i,]$column_import
    col_imp <- tbl_cols[which(grepl(col_imp, tbl_cols))][1]
    col_def <- tbl_config[i,]$column_defined

    # detect the split
    if(length(which(grepl(pattern, tbl_data[,col_imp] %>% unlist()))) == 0) {
      cli::cli_alert_info("... split {.emph {col_def}} with {.emph {pattern}} not detected")
      next()
    }

    prev_nrows <- tbl_data %>% nrow()

    if(col_def == 'protein'){
      cli::cli_alert_info("... created a {.emph protein_group} accounting")
      tbl_data <- tbl_data %>% dplyr::mutate(protein_cluster = paste0('pg', dplyr::row_number()))
    }

    tbl_data <- tbl_data %>% tidyr::separate_rows(dplyr::matches(paste0("^", col_imp, "$")), sep=pattern)
    cli::cli_alert_info("... split {.emph {col_def}} with {.emph {pattern}} resulting in {tbl_data %>% nrow() - prev_nrows} new rows")
  }
  return(tbl_data)
}

#' A helper function for importing peptide table data
#'
#' @param tbl_data a table of imported data
#' @param tbl_config a table of config values
#'
#' @return a tibble
#'
import_remove <- function(
    tbl_data = NULL,
    tbl_config = NULL
){

  # visible bindings
  pattern_remove <- NULL

  tbl_config <- tbl_config %>% dplyr::filter(!is.na(pattern_remove)) %>% as.data.frame()
  if(nrow(tbl_config) == 0) return(tbl_data)
  for(i in row.names(tbl_config)){
    pattern <- tbl_config[i,]$pattern_remove
    col_imp <- tbl_config[i,]$column_import
    col_def <- tbl_config[i,]$column_defined

    col <- which(grepl(col_imp, names(tbl_data)))
    w <- which(grepl(pattern, unlist(tbl_data[,col])))
    if(length(w) > 0) {tbl_data <- tbl_data[-w,]}
    cli::cli_alert_info("... removed {.emph {pattern}} from {.emph {col_def}} in {length(w)} rows")

  }
  return(tbl_data)
}

#' A helper function for importing peptide table data
#'
#' @param tbl_data a table of imported data
#' @param tbl_config a table of config values
#'
#' @return a tibble
#'
import_mbr <- function(
    tbl_data = NULL,
    tbl_config = NULL
){

  if(!'match_between_runs' %in% colnames(tbl_data)) {
    tbl_data$match_between_runs <- FALSE
  }

  if(nrow(tbl_config) == 0) return(tbl_data)

  cols_mbr <- tbl_config[which(!is.na(tbl_config$pattern_extract) & tbl_config$category == "impute"),]

  # sort out match_between_runs
  if(nrow(cols_mbr) > 0) {
    tbl_data$match_between_runs <- grepl(cols_mbr$pattern_extract[1], unlist(tbl_data$match_between_runs))
  }

  n_mbr <- length(which(tbl_data$match_between_runs == TRUE))
  if(n_mbr > 1) {
    cli::cli_alert_info("... match between runs accounting for {.emph {signif(n_mbr / nrow(tbl_data) * 100,3)}%} of data")
  } else {
    cli::cli_alert_info("... match between runs not found in data")
  }

  return(tbl_data)
}

#' A helper function for importing peptide table data
#'
#' @param tbl_data a table of imported data
#' @param tbl_config a table of config values
#'
#' @return a tibble
#'
import_validate <- function(
    tbl_data = NULL,
    tbl_config = NULL
){

  # visible bindings
  column_import <- NULL
  column_defined <- NULL

  cols_needed <- c("sample_file", "sample", "abundance_raw")
  cols_all <- set_vect(tbl_config)
  cols_quantitative <- set_vect(tbl_config, 'quantitative') %>% names()
  cols_annotations <- set_vect(tbl_config, 'annotation') %>% names()
  cols_accountings <- set_vect(tbl_config, 'accounting') %>% names()
  cols_identifiers <- set_vect(tbl_config, 'identifier') %>% names()
  cols_samples <- set_vect(tbl_config, 'sample') %>% names()
  # cols_impute <- set_vect(tbl_config, 'impute') %>% names()
  cols_mbr <- tbl_config[which(!is.na(tbl_config$pattern_extract) & tbl_config$category == "impute"),]
  # cols_filter <- set_vect(tbl_config, 'filter') %>% names()

  tbl_config <- tbl_config %>%
    dplyr::mutate(column_import = ifelse(is.na(column_import), column_defined, column_import))

  tbl_data_cols <- tbl_data %>% colnames()
  get_cols_all <- intersect(tbl_data_cols, names(cols_all))

  cols_notfound <- setdiff(cols_needed, get_cols_all)
  if(length(cols_notfound) > 0) {
    cli::cli_abort(".. import error, did not find columns for {.emph {cols_notfound}}")
  }

  get_cols_quantitative <- intersect(tbl_data_cols, cols_quantitative)
  get_cols_annotations <- intersect(tbl_data_cols, cols_annotations)
  get_cols_accountings <- intersect(tbl_data_cols, cols_accountings)
  get_cols_identifiers <- intersect(tbl_data_cols, cols_identifiers)
  get_cols_samples <- intersect(tbl_data_cols, cols_samples)

  cli::cli_div(theme = list(span.emph = list(color = "#ff4500")))
  if(length(get_cols_annotations) != length(cols_annotations)) {cli::cli_alert_warning("... not all `annotation` columns imported, missing {.emph {setdiff(cols_annotations, tbl_data_cols)}}")}
  if(length(get_cols_quantitative) != length(cols_quantitative)) {cli::cli_abort("... `quantitative` import error, did not find {.emph {cols_quantitative}}")}
  if(length(get_cols_identifiers) != length(cols_identifiers)) {cli::cli_abort("... `identifier` import error, missing {.emph {setdiff(cols_identifiers, tbl_data_cols)}}")}
  if(length(get_cols_samples) != length(cols_samples)) {cli::cli_abort("... `sample` import error, missing {.emph {setdiff(cols_samples, tbl_data_cols)}}")}

}
