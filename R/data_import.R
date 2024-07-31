#' A helper function for importing peptide table data
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
  num_psms <- NULL
  value <- NULL

  column_defined <- NULL
  pivot <- NULL
  column_import <- NULL
  accounting <- NULL
  account_str <- NULL
  protein_cluster <- NULL

  if(is.null(file_names)) {cli::cli_abort(c("x" = "No file paths indicated"))}
  if(is.null(platform)) {cli::cli_abort(c("x" = "No file platform indicated"))}
  if(is.null(analyte)) {cli::cli_abort(c("x" = "No file analyte indicated"))}

  if(!is.null(path)) {
    if(!file.exists(path)) {cli::cli_abort("Config file not found for {.emph {path}}")}
    path_to_config <- path
  } else {
    path_to_config <- glue::glue("{platform}_{analyte}") %>% path_to_package_data()
  }

  tbl_config <- path_to_config %>%
    readr::read_tsv(show_col_types = FALSE, progress = FALSE) %>%
    as.data.frame() %>%
    dplyr::filter(!is.na(column_defined))

  # does the config file have all the right columns
  tbl_config_need <- c("category","column_defined","column_import","pattern_extract","pattern_remove","pattern_split","pivot")
  tbl_config_have <- names(tbl_config)
  if(length(setdiff(tbl_config_need, tbl_config_have)) > 0) {
    cli::cli_div(theme = list(span.emph = list(color = "#ff4500")))
    cli::cli_abort(c("x" = "Columns missing in config fill: {.emph {setdiff(tbl_config_need, tbl_config_have)}}"))
  }

  ##############################################################################
  # order of operations
  # read_each >
  #        pivot_long > extract_pivot > pivot_wide >
  #        split > remove >
  #        extract >
  #        rename
  # merge_each
  # replicate sample accounting
  ##############################################################################

  # parse out tables for import manipulation
  tbl_pivot <- tbl_config %>% dplyr::filter(pivot == TRUE | is.na(column_import))
  pivot_cols <- tbl_pivot %>% dplyr::filter(!is.na(column_import)) %>% dplyr::select(column_import) %>% unlist() %>% as.character()
  pivot_str <- paste(pivot_cols, collapse = "|")

  # read in the config file for annotation columns
  cols_all <- set_vect(tbl_config)
  cols_all_str <- paste(cols_all[which(!is.na(cols_all))], collapse = "|")

  cli::cli_div(theme = list(span.emph = list(color = "blue")))
  cli::cli_alert_info("Importing {.emph {platform}}:")

  dat_out <- c()
  for(i in 1:length(file_names) ){
    file_name <- file_names[i]

    cli::cli_progress_step("... parsing {.emph {basename(file_name)}}")

    ############################################################################
    # read_data covers csv, tsv, xlsx and mzTab
    this_dat_obj <- file_name %>% read_data(platform, analyte)
    this_dat <- this_dat_obj$data %>% dplyr::select(dplyr::matches(cols_all_str, perl = TRUE))
    analyte <- this_dat_obj$analyte
    platform <- this_dat_obj$platform

    if(nrow(this_dat) == 0) {cli::cli_abort("... no data extracted, check config file")}

    if(pivot_str != "") {

      this_dat <- this_dat %>% dplyr::mutate(row_id = dplyr::row_number())

      this_pivot_not <- this_dat %>%
        dplyr::select(!dplyr::matches(pivot_str, perl = TRUE))

      for( pivot_col in pivot_cols ){

        pivot_col_name <- tbl_pivot %>%
          dplyr::filter(column_import == pivot_col) %>%
          dplyr::select('column_defined') %>%
          as.character()

        this_pivot <- this_dat %>%
          dplyr::select(dplyr::matches(paste0('row_id|', pivot_col), perl = TRUE)) %>%
          tidyr::pivot_longer(cols = dplyr::matches(pivot_col, perl = TRUE), names_to = 'accounting')

        # pivot extract
        tbl_pivot_extract <- tbl_pivot %>%
          dplyr::filter(is.na(column_import)) %>%
          dplyr::mutate(column_import = 'accounting')

        this_pivot <- this_pivot %>%
          import_extract(tbl_pivot_extract) %>%
          unique() %>%
          dplyr::select(!dplyr::matches('accounting')) %>%
          dplyr::rename(!!pivot_col_name := value)

        this_pivot_not <- this_pivot_not %>%
          dplyr::left_join(
            this_pivot,
            by = intersect(colnames(this_pivot_not), colnames(this_pivot))
          )
      }
      this_dat <- this_pivot_not
    }

    ############################################################################
    # split rows
    this_dat <- this_dat %>% import_split(tbl_config)

    ############################################################################
    # remove rows
    this_dat <- this_dat %>% import_remove(tbl_config)

    ############################################################################
    # extract values
    this_dat <- this_dat %>% import_extract(tbl_config, remove = TRUE)

    ############################################################################
    # rename columns
    this_dat <- this_dat %>% import_rename(tbl_config)

    ############################################################################
    # deal with protein_groups
    if(!'protein_cluster' %in% names(this_dat)) {
      cli::cli_abort("A `pattern_split` is required for identifier::protein in your config file")
    }
    u_prot_all <- this_dat$protein %>% unique() %>% length()
    this_dat <- this_dat %>%
      dplyr::group_by(protein_cluster) %>%
      dplyr::mutate(protein_group = paste(protein, collapse = "; ")) %>%
      dplyr::ungroup() %>%
      dplyr::select(!dplyr::matches('protein_cluster'))

    u_prot_red <- this_dat$protein %>% unique() %>% length()
    if(u_prot_red == u_prot_all){
      cli::cli_alert_info("... no homology detected")
    } else {
      cli::cli_alert_info("... homology collapsed {.emph {u_prot_all}} proteins into {.emph {u_prot_red}} proteins groups")
    }

    ############################################################################
    # munge MBR
    this_dat <- this_dat %>% import_mbr(tbl_config)

    ############################################################################
    # validate the imported data is good
    this_dat %>% import_validate(tbl_config)

    this_dat <- this_dat %>%
      dplyr::mutate(sample = gsub("\\-", "_", sample)) %>%
      dplyr::mutate(sample = gsub("\\/", "_", sample)) %>%
      dplyr::mutate(sample = tolower(sample)) %>%
      dplyr::mutate(sample = ifelse(grepl("^[0-9]", sample), paste0("s", sample), sample)) %>%
      dplyr::mutate(protein = trimws(protein)) %>%
      dplyr::mutate(import_file = basename(file_name)) %>%
      dplyr::mutate(sample_id = hash_vector(paste(import_file, sample_file, sample, i))) %>%
      dplyr::relocate(import_file, sample_file, sample_id) %>%
      as.data.frame() %>%
      dplyr::mutate(abundance_raw = ifelse(abundance_raw == 0, NA, abundance_raw))

    if(!'num_psms' %in% names(this_dat)) {this_dat <- this_dat %>% dplyr::mutate(num_psms = 1)}
    get_cols_num <- names(this_dat)[which(grepl('^num_|^abundance_', names(this_dat)))]
    for(get_col_num in get_cols_num){ this_dat[,get_col_num] <- as.numeric(this_dat[,get_col_num]) }

    # de-duplication !?
    # get_cols <- names(this_dat)[which(!grepl('abundance_raw|match_between_runs|num_', names(this_dat)))]
    # this_dat <- this_dat %>%
    #   dplyr::group_by(dplyr::across(get_cols)) %>%
    #   dplyr::slice_min(match_between_runs, with_ties = FALSE) %>%
    #   dplyr::summarise(num_psms = sum(num_psms),
    #                    match_between_runs = max(match_between_runs) == 1,
    #                    abundance_raw = sum(abundance_raw),
    #                    .groups = 'drop')

    dat_out <- dat_out %>% dplyr::bind_rows(this_dat)

    cli::cli_progress_done()
  }

  dat_out <- dat_out
  cols_annotations <- set_vect(tbl_config, 'annotation') %>% names()
  cols_identifiers <- set_vect(tbl_config, 'identifier') %>% names()

  cli::cli_progress_step("... tidying the data and finishing up")
  dl <- list(
    origin = platform,
    analyte = analyte,
    identifier = cols_identifiers,
    quantitative_source = 'raw',
    operations = list()
  ) %>% append(
    dat_out %>%
      codify(
        identifier = cols_identifiers,
        annotations = cols_annotations
      ))
  cli::cli_progress_done()

  return(dl)
}
