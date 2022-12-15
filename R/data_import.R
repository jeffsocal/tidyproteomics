#' A helper function for importing peptide table data from ProteomeDiscoverer
#'
#' @param file_names a character vector of file paths
#' @param platform a character string
#' @param analyte a character string
#' @param path a character string
#'
#' @return a tidyproteomics list data-object
#'
data_import <- function(
    file_names = NULL,
    platform = NULL,
    analyte = NULL,
    path = NULL
){

  # visible bindings
  protein <- NULL
  import_file <- NULL
  sample_file <- NULL
  sample_id <- NULL
  abundance_raw <- NULL
  match_between_runs <- NULL
  num_psms <- NULL

  if(is.null(file_names)) {cli::cli_abort(c("x" = "No file paths indicated"))}
  if(is.null(platform)) {cli::cli_abort(c("x" = "No file platform indicated"))}
  if(is.null(analyte)) {cli::cli_abort(c("x" = "No file analyte indicated"))}

  if(!is.null(path)) {
    if(!file.exists(path)) {cli::cli_abort("Config fine not found for {.emph {path}}")}
    path_to_config <- path
  } else {
    path_to_config <- glue::glue("{platform}_{analyte}") %>% path_to_package_data()
  }

  cols_config <- path_to_config %>%
    readr::read_tsv(show_col_types = FALSE, progress = FALSE) %>%
    as.data.frame()

  cols_needed <- c("sample_file", "sample", "abundance_raw")
  cols_extract <- cols_config[which(!is.na(cols_config$pattern_extract) & cols_config$category != "impute"),]
  cols_remove <- cols_config[which(!is.na(cols_config$pattern_remove)),]
  cols_mbr <- cols_config[which(!is.na(cols_config$pattern_extract) & cols_config$category == "impute"),]

  # read in the config file for annotation columns
  cols_all <- set_vect(cols_config)
  cols_quantitative <- set_vect(cols_config, 'quantitative')
  cols_annotations <- set_vect(cols_config, 'annotation')
  cols_accountings <- set_vect(cols_config, 'accounting')
  cols_identifiers <- set_vect(cols_config, 'identifier')
  cols_samples <- set_vect(cols_config, 'sample')
  cols_impute <- set_vect(cols_config, 'impute')
  cols_filter <- set_vect(cols_config, 'filter')

  dat_out <- c()
  for(i in 1:length(file_names) ){
    file_name <- file_names[i]

    this_dat <- file_name %>% read_data(show_col_types = FALSE, progress = FALSE)

    cli::cli_progress_step("Importing {.emph {basename(file_name)}}")

    this_dat_cols <- this_dat %>% colnames()
    get_cols_all <- match_vect(this_dat_cols, cols_all)

    cols_notfound <- setdiff(cols_needed, names(get_cols_all))
    if(length(cols_notfound) > 0) {
      cli::cli_abort(".. import error, did not find columns for {.emph {cols_notfound}}")
    }

    get_cols_quantitative <- match_vect(this_dat_cols, cols_quantitative)
    get_cols_annotations <- match_vect(this_dat_cols, cols_annotations)
    get_cols_accountings <- match_vect(this_dat_cols, cols_accountings)
    get_cols_identifiers <- match_vect(this_dat_cols, cols_identifiers)
    get_cols_samples <- match_vect(this_dat_cols, cols_samples)
    get_cols_impute <- match_vect(this_dat_cols, cols_impute)
    get_cols_filter <- match_vect(this_dat_cols, cols_filter)

    if(length(get_cols_annotations) != length(cols_annotations)) {cli::cli_alert_warning(".. not all `annotation` columns imported, only found {.emph {get_cols_annotations}}")}
    if(length(get_cols_quantitative) != length(cols_quantitative)) {cli::cli_abort(".. `quantitative` import error, did not find {.emph {cols_quantitative}}")}
    if(length(get_cols_identifiers) != length(cols_identifiers)) {cli::cli_abort(".. `identifier` import error, only found {.emph {get_cols_identifiers}}")}
    if(length(get_cols_samples) != length(cols_samples)) {cli::cli_abort(".. `sample` import error, only found {.emph {get_cols_samples}}")}

    this_dat <- this_dat %>%
      dplyr::select(dplyr::all_of(c(get_cols_identifiers,
                                    get_cols_all))) %>%
      tidyr::separate_rows(protein, sep = "\\;\\s*") %>%
      dplyr::mutate(sample = gsub("\\-", "_", sample)) %>%
      dplyr::mutate(sample = gsub("\\/", "_", sample)) %>%
      dplyr::mutate(sample = ifelse(grepl("^[0-9]", sample), paste0("s", sample), sample)) %>%
      dplyr::mutate(protein = trimws(protein)) %>%
      dplyr::mutate(import_file = basename(file_name)) %>%
      dplyr::mutate(sample_id = hash_vector(paste(import_file, sample_file, i))) %>%
      dplyr::relocate(import_file, sample_file, sample_id) %>%
      as.data.frame()

    # subset the data table based on patterns to remove
    if(nrow(cols_remove) > 0) {
      for(i in 1:nrow(cols_remove)) {
        w <- which(grepl(cols_remove$pattern_remove[i], this_dat[,cols_remove$column_defined[i]]))
        if(length(w) > 0) {this_dat <- this_dat[-w,]}
      }
    }

    # pull values out of the data based on a pattern
    if(nrow(cols_extract) > 0) {
      for(i in 1:nrow(cols_extract)) {
        this_dat[,cols_extract$column_defined[i]] <- stringr::str_extract(this_dat[,cols_extract$column_defined[i]], cols_extract$pattern_extract[i])
      }
    }

    # remove filtering columns
    if(length(cols_filter) > 0) {this_dat[,names(cols_filter)] <- NULL}

    # sort out match_between_runs
    if(nrow(cols_mbr) > 0) {
      this_dat$match_between_runs <- grepl(cols_mbr$pattern_extract[1], this_dat$match_between_runs)
    } else {
      this_dat$match_between_runs <- FALSE
    }

    get_cols <- names(this_dat)[which(!grepl('abundance_raw|match_between_runs|num_', names(this_dat)))]

    if(!'num_psms' %in% names(this_dat)) {this_dat <- this_dat %>% dplyr::mutate(num_psms = 1)}

    this_dat <- this_dat %>%
      dplyr::mutate(abundance_raw = as.numeric(abundance_raw)) %>%
      dplyr::group_by(dplyr::across(get_cols)) %>%
      dplyr::slice_max(match_between_runs) %>%
      dplyr::summarise(num_psms = sum(num_psms),
                       match_between_runs = max(match_between_runs) == 1,
                       abundance_raw = sum(abundance_raw),
                       .groups = 'drop') %>%
      dplyr::mutate(abundance_raw = ifelse(abundance_raw == 0, NA, abundance_raw))

    dat_out <- dat_out %>% dplyr::bind_rows(this_dat)

    cli::cli_progress_done()
  }

  dl <- list(
    origin = platform,
    analyte = analyte,
    identifier = names(cols_identifiers),
    quantitative_source = 'raw',
    operations = list()
  ) %>% append(
    dat_out %>%
      codify(
        identifier = names(cols_identifiers),
        annotations = names(cols_annotations)
      ))

  return(dl)
}
