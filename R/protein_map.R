#' Align a peptide data to protein sequences for visualization
#'
#' @param data a tidyproteomics data-object, specifically of peptide origin
#' @param fasta_path a character string representing the path to a fasta file
#'
#' @return a list of protein mappings
#' @export
#'
#' @examples
#' library(dplyr, warn.conflicts = FALSE)
#' library(tidyproteomics)
#'
#' hela_protein_map <- hela_peptides %>%
#'    protein_map(fasta = path_to_package_data('fasta'))
#'
protein_map <- function(
    data = NULL,
    fasta_path = NULL
){

  # visible bindings
  identifier <- NULL
  protein <- NULL
  peptide <- NULL
  modifications <- NULL
  mod_map <- NULL
  position <- NULL
  position_start <- NULL
  probability <- NULL
  modification <- NULL
  aa <- NULL
  abundance <- NULL
  frequency <- NULL
  position_end <- NULL

  check_data(data)
  l_fasta <- rfasta::parse(fasta_path)
  data_quant <- tidyproteomics::extract(data, data$quantitative_source) %>%
    tidyr::separate(identifier, into = c('protein', 'peptide', 'modifications'),
                    sep = "\\s\\|\\s")

  vrange <- function(x, y){ return(x:y)}
  min_abundance <- min(data_quant$abundance, na.rm = T)

  # collect all the proteins
  proteins <- data_quant %>%
    dplyr::select(protein) %>%
    unique() %>%
    unlist() %>%
    as.character()

  cli::cli_process_start("Sequence mapping")
  cli::cli_progress_bar("  Sequence mapping", total  = proteins %>% length())

  l_out <- list()
  for( protein in proteins ){

    cli::cli_progress_update()

    # pull out the protein of interest
    l_fasta_this <- l_fasta[[protein]]
    if(length(l_fasta_this) == 0) {next()}

    this_protein <- l_fasta_this$sequence
    this_peptide <- data_quant %>% dplyr::filter(protein == protein)

    peptide_map <- align_peptide(this_peptide$peptide %>% unique(), this_protein) %>%
      dplyr::inner_join(this_peptide, by = 'peptide')

    modification_map <- peptide_map %>%
      dplyr::mutate(mod_map = purrr::map2(peptide, modifications, align_modification)) %>%
      tidyr::unnest(mod_map) %>%
      dplyr::mutate(position = position + position_start - 1) %>%
      dplyr::mutate(probability = ifelse(is.na(probability), 1, probability)) %>%
      tidyr::separate_rows(modification, sep="; ") %>%
      dplyr::group_by(aa, position, peptide, modifications,
                      modification, sample, replicate) %>%
      dplyr::summarise(
        abundance = sum(abundance, na.rm = T),
        .groups = "drop"
      ) %>%
      dplyr::group_by(aa, position, sample, replicate) %>%
      dplyr::mutate(frequency = abundance / sum(abundance, na.rm = T)) %>%
      dplyr::ungroup() %>%
      dplyr::select(aa, position, peptide, modifications, modification,
                    sample, replicate, abundance, frequency)

    # build a map of the protein for each sample-replicate
    protein_map <- tibble::tibble(
      aa = stringr::str_extract_all(this_protein, "[A-Z]{1}")[[1]]) %>%
      dplyr::mutate(position = dplyr::row_number()) %>%
      base::merge(this_peptide %>% dplyr::select(sample, replicate) %>% unique())

    peptide_map_position <- peptide_map %>%
      dplyr::mutate(position = purrr::map2(position_start, position_end, vrange)) %>%
      tidyr::unnest(position)

    protein_map <- protein_map %>%
      dplyr::full_join(peptide_map_position, by = c('position', 'sample','replicate')) %>%
      dplyr::group_by(aa, position, sample, replicate) %>%
      dplyr::summarise(
        proteins = length(position[abundance > 0 & !is.na(abundance)]),
        abundance = sum(abundance, na.rm = T),
        .groups = 'drop'
      ) %>%
      dplyr::arrange(position)

    protein_map <- protein_map %>% dplyr::group_by(aa, position, sample) %>%
      dplyr::summarise(abundance = stats::median(abundance, na.rm=T), .groups = 'drop') %>%
      dplyr::mutate(protein = protein) %>% dplyr::relocate(protein) %>%
      dplyr::arrange(position)

    peptide_map <- peptide_map %>% dplyr::group_by(peptide, position_start, position_end,
                                            modifications, sample) %>%
      dplyr::summarise(abundance = stats::median(abundance, na.rm=T), .groups = 'drop') %>%
      dplyr::mutate(protein = protein) %>% dplyr::relocate(protein) %>%
      dplyr::arrange(position_start)

    modification_map <- modification_map %>% dplyr::group_by(aa, position, modification, sample) %>%
      dplyr::summarise(abundance = stats::median(abundance, na.rm=T), .groups = 'drop') %>%
      dplyr::mutate(protein = protein) %>% dplyr::relocate(protein) %>%
      dplyr::arrange(position) %>%
      dplyr::filter(modification != "")


    l_out[[protein]] <- list(
      protein = protein,
      sequence = this_protein,
      coverage = protein_map %>%
        dplyr::group_by(sample) %>%
        dplyr::summarise(coverage = position %>% unique() %>% length() / stringr::str_length(this_protein), .groups='drop'),
      residues =  protein_map,
      peptides = peptide_map,
      modifications = modification_map
    )

  }

  cli::cli_process_done

  l_out <- list(
    "origin" = data$origin,
    "analyte" = "sequences",
    "quantitative" = l_out
  )

  class(l_out) <- 'tidyproteomics'

  return(l_out)
}
