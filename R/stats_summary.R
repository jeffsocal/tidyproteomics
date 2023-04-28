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
  identifier <- data$identifier

  cli::cli_progress_cleanup()
  cli::cli_progress_bar(type = 'tasks')
  if(data$analyte == 'peptides'){
    cli::cli_progress_step(" ... calculating peptide -to- protein stats")
    data <- data %>% collapse(.verbose = FALSE)
  }

  tbl_long <- data %>% meld(single_quant_source = TRUE) %>%
    munge_identifier('combine', identifiers = identifier) %>%
    dplyr::mutate(identifiers = stringr::str_count(identifier, ";") + 1)

  tbl_stats <- tbl_long %>%
    dplyr::group_by_at(group_by) %>%
    dplyr::summarise(
      identifier_groups = sum(identifiers),
      identifiers = identifier %>% unique() %>% length(),
      peptides = ceiling(num_peptides %>% sum(na.rm = T) / length(unique(sample_id))),
      peptides_unique = ceiling(num_unique_peptides %>% sum(na.rm = T) / length(unique(sample_id))),
      quantifiable = signif(num_peptides[!is.na(abundance)] %>% sum(na.rm=T) / length(unique(sample_id)) / peptides * 100, 3),
      .groups = 'drop'
    )

  tbl_cvs <- tbl_long %>%
    dplyr::group_by_at(c(group_by, 'identifier')) %>%
    dplyr::summarise(
      CVs = signif(sd(abundance, na.rm = TRUE)/mean(abundance, na.rm = TRUE),2),
      .groups = 'drop'
    ) %>%
    dplyr::group_by_at(group_by) %>%
    dplyr::summarise(
      CVs = median(CVs, na.rm = TRUE),
      .groups = 'drop'
    ) %>%
    dplyr::filter(!is.na(CVs))

  if(nrow(tbl_cvs) != 0) { tbl_stats <- bind_cols(tbl_stats, tbl_cvs) }

  w_cols <- which(grepl('^identifier[s_]*', colnames(tbl_stats)))
  colnames(tbl_stats)[w_cols] <- sub('identifier', identifier, colnames(tbl_stats)[w_cols])

  cli::cli_process_done()
  cli::cli_progress_cleanup()

  return(tbl_stats)
}


