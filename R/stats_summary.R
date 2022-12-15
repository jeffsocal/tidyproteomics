#' Summarize the protein accounting
#'
#' @description
#' `stats_summary()` is an analysis function that computes the protein summary
#' statistics for a given tidyproteomics data object.
#'
#' @param data tidyproteomics data object
#' @param group_by what to summarize
#'
#' @return a tibble
#'
stats_summary <- function(
    data,
    group_by = c("global", "sample", "replicate", "experiment")
){

  # visible bindings
  sample_id <- NULL
  files <- NULL
  peptide <- NULL
  abundance <- NULL
  matches <- NULL
  identifier <- NULL
  sd <- NULL
  median <- NULL
  cv <- NULL
  num_peptides <- NULL
  num_unique_peptides <- NULL
  peptides <- NULL
  quantifiable <- NULL
  protein <- NULL
  abundance_raw <- NULL
  protein_group <- NULL
  proteins <- NULL

  group_by <- rlang::arg_match(group_by)
  if(group_by == 'experiment') { group_by <- 'sample_id' }
  if(group_by == 'global') { group_by <- NULL }
  check_data(data)

  cli::cli_progress_cleanup()
  cli::cli_progress_bar(type = 'tasks')
  if(data$analyte == 'peptides'){
    cli::cli_progress_step(" ... calculating peptide -to- protein stats")
    data <- data %>% collapse(.verbose = FALSE)
  }

  data_quant <- data$quantitative %>%
    dplyr::rename(abundance = abundance_raw, identifier = protein)
  data_exps <- data$experiments
  data_account <- data$accounting %>% dplyr::full_join(data_exps, by = c('sample_id')) %>%
    dplyr::rename(identifier = protein)

  tb_temp_e <- data_exps %>%
    dplyr::group_by(dplyr::across(c(group_by))) %>%
    dplyr::summarise(files = dplyr::n(), .groups = 'drop')

  tb_temp_qnt <- data_account %>%
    dplyr::inner_join(data_quant,
                      by = c("identifier", "sample_id", "sample", "replicate")) %>%
    dplyr::mutate(proteins = stringr::str_count(protein_group, ";") + 1) %>%
    dplyr::group_by(dplyr::across(c(group_by))) %>%
    dplyr::summarise(
      protein_groups = protein_group %>% unique() %>% length(),
      proteins = identifier %>% unique() %>% length(),
      peptides = ceiling(num_peptides %>% sum(na.rm = T) / length(unique(sample_id))),
      peptides_unique = ceiling(num_unique_peptides %>% sum(na.rm = T) / length(unique(sample_id))),
      quantifiable = signif(num_peptides[!is.na(abundance)] %>% sum(na.rm=T) / length(unique(sample_id)) / peptides * 100, 3),
      .groups = 'drop'
    )
  cli::cli_process_done()
  cli::cli_progress_cleanup()

  if(!is.null(group_by) && group_by == 'sample_id'){
    tb_temp_qnt <- data_exps %>%
      dplyr::inner_join(tb_temp_qnt, by = 'sample_id')

    return(tb_temp_qnt)
  }

  tb_temp_a <- data_quant %>%
    dplyr::group_by(dplyr::across(c(group_by, identifier))) %>%
    dplyr::summarise(
      cv = sd(abundance, na.rm=T) / mean(abundance, na.rm=T),
      .groups = 'drop'
    ) %>%
    dplyr::group_by(dplyr::across(c(group_by))) %>%
    dplyr::summarise(
      CVs = median(cv, na.rm=T),
      .groups = 'drop'
    )

  if(!is.null(group_by)) {
    tb_temp_qnt <- tb_temp_qnt %>%
      dplyr::full_join(tb_temp_a, by = group_by) %>%
      dplyr::full_join(tb_temp_e, by = group_by) %>%
      dplyr::relocate(files)
  } else {
    tb_temp_qnt <- tb_temp_qnt %>%
      dplyr::bind_cols(tb_temp_a) %>%
      dplyr::bind_cols(tb_temp_e) %>%
      dplyr::relocate(files)
  }

  return(tb_temp_qnt)
}


