#' Align a peptide data to protein sequences for visualization
#'
#' @param mapped_data a tidyproteomics data-object, specifically of peptide origin
#' @param protein a character string
#' @param row_length a numeric
#' @param samples a character string
#' @param modifications a character string
#'
#' @return a plot munged list of protein mappings
#'
protein_map_munge <- function(
    mapped_data = NULL,
    protein = NULL,
    row_length = 50,
    samples = NULL,
    modifications = NULL
){

  # visible bindings
  coverage <- NULL
  abundance <- NULL
  modification <- NULL
  position <- NULL
  position_start <- NULL
  position_end <- NULL
  row_a <- NULL
  row_b <- NULL
  position_end_a <- NULL
  position_start_a <- NULL
  peptide <- NULL
  row_n <- NULL
  mod_map <- NULL
  abundance_residue <- NULL

  check_data(mapped_data)

  protein <- rlang::arg_match(protein, names(mapped_data$quantitative))
  mapped_data <- mapped_data$quantitative[[protein]]
  if(mode(row_length) != 'numeric') {cli::cli_abort(c("i" = "`row_length` must be a numeric, not a {mode(row_length)}"))}
  row_length <- floor(row_length)
  if(row_length < 10 || row_length > 500) {cli::cli_abort(c("i" = "`row_length` must be between 10 and 500, {.emph {row_length}} is not accepted"))}
  if(!is.null(samples)) {samples <- rlang::arg_match(samples, mapped_data$covereage$sample, multiple = T)}
  if(!is.null(modifications)) {modifications <- rlang::arg_match(modifications, unique(mapped_data$mod_map$modification), multiple = T)}

  coverag_map <- mapped_data$coverage %>% dplyr::mutate(coverage = paste0(signif(coverage * 100, 3), "%"))

  residue_map <- mapped_data$residues %>% dplyr::mutate(abundance = ifelse(is.na(abundance), 0, abundance)) %>%
    dplyr::full_join(coverag_map, by = 'sample') %>% tidyr::unite(sample, sample, coverage, sep = " ")

  peptide_map <- mapped_data$peptides %>%
    dplyr::full_join(coverag_map, by = 'sample') %>% tidyr::unite(sample, sample, coverage, sep = " ")

  modification_map <- mapped_data$modifications %>%
    dplyr::full_join(coverag_map, by = 'sample') %>% tidyr::unite(sample, sample, coverage, sep = " ")

  n_mods <- modification_map %>% dplyr::select(modification) %>% unique() %>% base::nrow()

  munged_residues <- residue_map %>%
    dplyr::group_by(sample) %>%
    dplyr::mutate(
      col = rep(1:row_length,ceiling(max(position)/row_length))[1:max(position)],
      row = sort(rep(0:ceiling(max(position)/row_length),row_length), decreasing = F)[0:max(position)] * row_length + 1
    ) %>%
    dplyr::ungroup()

  munged_peptides <- peptide_map %>%
    dplyr::mutate(row_a = ceiling(position_start / row_length),
                  row_b = ceiling(position_end / row_length),
                  position_start_a = position_start - (row_a-1) * row_length,
                  position_end_a = position_end - (row_b-1) * row_length) %>%
    dplyr::mutate(position_end_a = ifelse(position_end_a < position_start_a, row_length, position_end_a))

  munged_peptides <- munged_peptides%>%
    dplyr::bind_rows(
      munged_peptides %>%
        dplyr::filter(row_a != row_b)  %>%
        dplyr::mutate(position_start_a = 1,
                      position_end_a = position_end - (row_b-1) * row_length,
                      row_a = row_b)
    )

  # segment mapping
  munged_peptides_byrow <- munged_peptides %>%
    dplyr::arrange(row_a, position_start_a) %>%
    dplyr::select(peptide) %>% unique() %>%
    dplyr::mutate(row_n = dplyr::row_number())

  munged_peptides_bysegment <- munged_peptides %>%
    dplyr::filter(!is.na(abundance) & abundance > 0) %>%
    dplyr::inner_join(munged_peptides_byrow, by = 'peptide') %>%
    dplyr::group_by(peptide, modifications, sample, position_start_a, position_end_a, row_a, row_n) %>%
    dplyr::summarise(abundance = sum(abundance, na.rm=T), .groups = "drop") %>%
    dplyr::group_by(sample) %>%
    dplyr::mutate(row = (row_a - 1) * row_length +  row_length * (row_n/max(munged_peptides_byrow$row_n))) %>%
    dplyr::ungroup()

  munged_modifications_bysegment <- munged_peptides_bysegment %>%
    dplyr::mutate(mod_map = purrr::map2(peptide, modifications, align_modification)) %>%
    tidyr::unnest(mod_map) %>%
    dplyr::filter(position >= position_start_a,
                  position <= position_end_a)

  munged_modifications <- munged_residues %>%
    dplyr::rename(abundance_residue = abundance) %>%
    dplyr::inner_join(modification_map, by = c('protein','aa','position','sample')) %>%
    dplyr::arrange(modification) %>%
    dplyr::group_by(sample, position) %>%
    dplyr::mutate(row_d = dplyr::row_number() + (row_length*.6) / n_mods,
                  frequency = abundance / abundance_residue) %>%
    dplyr::ungroup()

  if(!is.null(modifications)) {
    munged_modifications <- munged_modifications %>% dplyr::filter(modification %in% modifications)
  }

  if( !is.null(samples) ){

    msr <- residue_map$sample %>% unique()
    msr <- msr[which(grepl(samples, msr))]

    munged_residues <- munged_residues %>% dplyr::filter(sample %in% msr)
    munged_modifications <- munged_modifications %>% dplyr::filter(sample %in% msr)

    munged_peptides_bysegment <- munged_peptides_bysegment %>% dplyr::filter(sample %in% msr)
    munged_modifications_bysegment <- munged_modifications_bysegment %>% dplyr::filter(sample %in% msr)
  }

  l_out <- list(
    munged_residues = munged_residues %>% dplyr::arrange(position),
    munged_peptides = munged_peptides %>% dplyr::arrange(position_start),
    munged_modifications = munged_modifications %>% dplyr::arrange(position),
    segments = list(
      peptides = munged_peptides_bysegment,
      modifications = munged_modifications_bysegment
    )
  )

  return(l_out)
}
