#' Align a peptide sequence to a protein sequence
#'
#' @param peptide a character string representing a peptide sequence
#' @param protein a character string representing a protein sequence
#'
#' @return a tidyproteomics data-object
#'
align_peptide <- function(
    peptide = NULL,
    protein = NULL
){

  if(mode(peptide) != "character") {cli::cli_abort(c("x" = "peptide is `{mode(peptide)}`, should be a character"))}
  if(mode(protein) != "character") {cli::cli_abort(c("x" = "protein is `{mode(protein)}`, should be a character"))}

  if(!is.vector(peptide)) { peptide <- c(peptide) }

  l_m <- list()
  for(i in 1:length(peptide)) {
    mp <- Biostrings::matchPattern(peptide[i], protein)

    s <- mp@ranges@start
    w <- mp@ranges@width

    l_m[[i]] <- tibble::tibble(
      'peptide' = peptide[i],
      'position_start' = s,
      'position_end' = s + w -1
    )
  }
  l_m %>% dplyr::bind_rows()
}

#' Align a modification to a peptide sequence
#'
#' @param peptide a character string representing a peptide sequence
#' @param modification a character string representing a modification and location probability
#'
#' @return a tidyproteomics data-object
#'

#'
align_modification <- function(
    peptide = NULL,
    modification = NULL
){

  # visible bindings
  aa <- NULL
  position <- NULL
  probability <- NULL

  pep_map <- tibble::tibble(aa = stringr::str_extract_all(peptide, "[A-Z]{1}")[[1]]) %>%
    dplyr::mutate(position = dplyr::row_number())

  # if there is no modification return as is
  if(is.na(modification) || !grepl("\\[", modification)){
    return(pep_map %>% dplyr::mutate(modification = '', probability = NA))
  }

  l_mods <- stringr::str_split(modification, "\\]\\;")[[1]]

  mod_dt <- list()
  mod_dn <- 0

  for(mod in l_mods){
    mod <- trimws(mod)
    mod_n <- as.integer(stringr::str_extract(mod, '^\\d+'))
    mod_t <- sub("\\s\\[.+", "", sub(paste0(mod_n, 'x'), '', mod))
    mod_as <- stringr::str_split(sub("\\]", "", sub(".+\\s\\[", "", mod)),";")
    mod_as <- mod_as[[1]]

    mod_p <- 1

    if(mod_as %>% length() != mod_n) { mod_p = 1/mod_n }

    for( mod_a in mod_as){

      mod_a_n <- stringr::str_extract(mod_a, '\\d+')
      mod_a_s <- stringr::str_split(gsub("\\d+", "", mod_a), '/')

      mod_p = mod_p/length(mod_a_s[[1]])

      for( mod_a_a in mod_a_s[[1]] ){

        mod_dn <- mod_dn + 1;

        mod_pt <- mod_p
        if(is.na(mod_a_n)) { mod_pt <- mod_p / pep_map %>% dplyr::filter(aa == mod_a_a) %>% nrow() }

        mod_dt[[mod_dn]] <- tibble::tibble(
          aa = trimws(mod_a_a),
          position = mod_a_n %>% as.integer(),
          modification = mod_t,
          probability = mod_pt
        )

      }
    }

    mod_map <- mod_dt %>% dplyr::bind_rows()

  }

  mod_map_ind <- mod_map %>% dplyr::filter(!is.na(position))
  mod_map_unk <- mod_map %>% dplyr::filter(is.na(position))

  pep_map_ind <- pep_map_unk <- pep_map

  if(mod_map_ind %>% nrow() != 0){
    pep_map_ind <- pep_map_ind %>%
      dplyr::full_join(mod_map_ind, by = c('aa','position'))
  }

  if(mod_map_unk %>% nrow() != 0){
    pep_map_unk <- pep_map_unk %>%
      dplyr::full_join(mod_map_unk %>% dplyr::select(-position), by = c('aa'))
  }

  pep_map <- pep_map_ind %>%
    dplyr::bind_rows(pep_map_unk) %>%
    dplyr::filter(!is.na(modification)) %>%
    dplyr::group_by(aa, position) %>%
    dplyr::summarise(
      modification = paste(modification, collapse = '; '),
      probability = mean(probability, na.rm = T),
      .groups = 'drop'
    ) %>%
    dplyr::mutate(modification = gsub("\\;\\s$", "",  modification),) %>%
    dplyr::arrange(position) %>%
    dplyr::filter(!is.na(position))

  return(pep_map)

}
