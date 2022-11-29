#' A helper function for importing peptide table data from ProteomeDiscoverer
#'
#' @param file_names a character vector of file paths
#'
#' @return a tidyproteomics list data-object
#'
import_pd_peps <- function(
    file_names = NULL
){

  # visible bindings
  .data <- NULL
  abundance_raw <- NULL
  annotation <- NULL
  sample_file <- NULL
  sample_file_mbr <- NULL
  match_between_runs <- NULL
  peptide <- NULL
  protein <- NULL
  sample <- NULL
  term <- NULL
  import_file <- NULL

  if(is.null(file_names)) {
    cli::cli_abort(c("x" = "No file paths indicated"))
  }

  cols_config <- path_to_package_data("proteomediscoverer_peptide_headers") %>%
    readr::read_tsv(col_names = c('category', 'column_defined', 'column_import'), show_col_types = FALSE, progress = FALSE) %>%
    as.data.frame()

  # read in the config file for annotation columns
  cols_quantitative <- set_vect(cols_config, 'quantitative')
  cols_annotations <- set_vect(cols_config, 'annotation')
  cols_accountings <- set_vect(cols_config, 'accounting')
  cols_identifiers <- set_vect(cols_config, 'identifier')
  cols_samples <- set_vect(cols_config, 'sample')
  cols_impute <- set_vect(cols_config, 'impute')

  dat_out <- c()

  cli::cli_progress_bar(type = 'tasks')
  for(i in 1:length(file_names) ){
    file_name <- file_names[i]

    this_dat <- file_name %>% readxl::read_xlsx()

    cli::cli_progress_step(paste("Importing {.emph {basename(file_name)}}"))

    tryCatch({

      this_dat_cols <- this_dat %>% colnames()
      get_cols_quantitative <- this_dat_cols[which(grepl(cols_quantitative, this_dat_cols))]
      get_cols_annotations <- match_vect(this_dat_cols, cols_annotations)
      get_cols_accountings <- match_vect(this_dat_cols, cols_accountings)
      get_cols_identifiers <- match_vect(this_dat_cols, cols_identifiers)
      get_cols_samples <- match_vect(this_dat_cols, cols_samples)
      get_cols_impute <- this_dat_cols[which(grepl(cols_impute, this_dat_cols))]

      this_dat <- this_dat %>%
        dplyr::select(dplyr::all_of(c(get_cols_identifiers,
                                      get_cols_samples,
                                      get_cols_quantitative,
                                      get_cols_impute,
                                      get_cols_accountings,
                                      get_cols_annotations))) %>%
        tidyr::separate_rows(protein, sep = "\\;\\s*") %>%
        dplyr::mutate(peptide = gsub("(\\[.+?\\]|\\.)", "", peptide)) %>%
        dplyr::mutate(protein = sub("^..\\|", "", protein)) %>%
        dplyr::mutate(protein = sub("\\|.+$", "", protein))  %>%
        tidyr::pivot_longer(dplyr::matches(cols_quantitative[[1]]),
                            names_to = 'sample', values_to = 'abundance_raw') %>%
        dplyr::mutate(sample_file = sub(cols_quantitative[[1]], "", sample) %>% tolower()) %>%
        dplyr::mutate(sample_file = sub("\\:.*", "", sample_file) %>% trimws()) %>%
        dplyr::mutate(import_file = basename(file_name)) %>%
        dplyr::mutate(sample_id = hash_vector(paste(import_file, sample_file, i))) %>%
        dplyr::mutate(sample = sub(".*(Sample|Control),\\s", "", sample) %>% tolower()) %>%
        dplyr::mutate(sample = sub("\\-", "_", sample)) %>%
        dplyr::mutate(sample = gsub("\\/", "_", sample)) %>%
        dplyr::mutate(abundance_raw = ifelse(abundance_raw == 0, NA, abundance_raw))

      if(length(get_cols_impute) > 0){
        this_dat <- this_dat %>%
          tidyr::pivot_longer(dplyr::matches(cols_impute[[1]]),
                              names_to = 'sample_file_mbr', values_to = 'match_between_runs') %>%
          dplyr::mutate(sample_file_mbr = stringr::str_extract(sample_file_mbr, "F[0-9]*\\:\\s[0-9]{0,3}[CN]*") %>% tolower()) %>%
          dplyr::mutate(sample_file_mbr = sub("\\:", "", sample_file_mbr) %>% trimws()) %>%
          dplyr::mutate(sample_file_mbr = sub("\\s+", "_", sample_file_mbr)) %>%
          dplyr::mutate(sample_file_mbr = sub("\\s/", "_", sample_file_mbr)) %>%
          dplyr::filter(sample_file_mbr == sample_file) %>%
          dplyr::select(!'sample_file_mbr') %>%
          dplyr::mutate(match_between_runs = ifelse(is.na(match_between_runs), "Not Found", match_between_runs)) %>%
          dplyr::mutate(match_between_runs = match_between_runs != 'High')
      } else {
        this_dat <- this_dat %>% dplyr::mutate(match_between_runs = FALSE)
      }

      # remove duplicates
      this_dat <- this_dat %>% unique()

      # filter out MBR with 0 quant value as this should indicate no-evidence
      w_mbr <- which(this_dat$match_between_runs == TRUE & this_dat$abundance_raw == 0)
      if(length(w_mbr) > 0) {this_dat <- this_dat[-w_mbr,]}

    }, error = function(err) {
      cli::cli_abort(c("x" = as.character(as.vector(err))))
    })

    if( length(get_cols_annotations) != length(cols_annotations) )
      cli::cli_alert_warning(".. not all expected annotation columns present in import file")

    dat_out <- dat_out %>% dplyr::bind_rows(this_dat)

  }
  cli::cli_progress_done()

  dl <- list(
    origin = 'ProteomeDiscoverer',
    analyte = 'peptides',
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
