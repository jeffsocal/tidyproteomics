#' A helper function for importing peptide table data
#'
#' @param path a character string
#' @param analyte a character string
#'
#' @return a tidyproteomics list data-object
#'
read_mzTab <- function(
    path = NULL,
    analyte = c("peptides","proteins")
){

  # visible bindings
  tbl_mtd <- NULL
  var <- NULL
  val <- NULL
  org <- NULL
  tbl_files <- NULL
  modifications <- NULL
  spectra_ref <- NULL
  charge <- NULL
  calc_mass_to_charge <- NULL
  num_psms <- NULL
  mass_to_charge <- NULL
  accession <- NULL

  rfn <- readr::read_lines(path)

  lcom <- rfn[which(grepl("^COM", rfn))]
  lmtd <- rfn[which(grepl("^MTD", rfn))]

  str_mtd <- function(x){
    x <- stringr::str_split(x, pattern ="\t") %>% unlist()
    x <- tibble::tibble(var = x[2], val = x[3])
    return(x)
  }

  vec_tbl <- function(x){
    x <- stringr::str_split(x, pattern ="\t") %>% unlist()
    names(x) <- paste0("V", 1:length(x))
    return(x)
  }

  tbl_names <- function(tbl, tbl_files) {
    col_names <- as.character(unlist(tbl[1,]))
    for(ii in 1:nrow(tbl_files)){
      assay_id <- tbl_files$assay_id[ii]
      sample_file <- tbl_files$sample_file[ii]
      sample_name <- tbl_files$sample_name[ii]
      col_names <- stringr::str_replace_all(col_names,
                                            glue::glue("\\[{assay_id}\\]"),
                                            glue::glue("[{sample_file}; {sample_name}]"))
    }

    colnames(tbl) <- col_names
    tbl <- tbl[-1,]
    return(tbl)
  }

  tbl_mtd <- lmtd %>% lapply(str_mtd) %>% dplyr::bind_rows()

  is_quant <- tbl_mtd %>%
    dplyr::filter(var == 'mzTab-type') %>%
    dplyr::select(val) %>% unlist() %>% as.character()

  if(is_quant != "Quantification") {
    cli::cli_div(theme = list(span.emph = list(color = "#ff4500")))
    cli::cli_alert_warning("mzTab file is not a type `Quantification`")
    cli::cli_abort("... is type {.emph {is_quant}}")
  }

  platform <- lcom[which(grepl('Quantification report',lcom))] %>%
    stringr::str_extract("(?<=report).+") %>%
    stringr::str_replace_all('^\\"\\s+|\\.$', "")

  # MODIFICATIONS:
  # generate a table of modifications for substitution into the PSM table
  tbl_mods <- tbl_mtd %>% dplyr::filter(grepl("\\_mod\\[", var)) %>%
    tidyr::separate(val, sep="\\,\\s*", into = c('org', 'id', 'str', 'extra')) %>%
    dplyr::mutate(str = ifelse(org == '[CHEMMOD', sub("CHEMMOD\\:", "", id), str)) %>%
    dplyr::select(id, str)

  # DATA FILES
  # generate a table of data files used in the analysis
  tbl_files <- tbl_mtd %>%
    dplyr::filter(grepl("^ms_run", var)) %>%
    dplyr::mutate(var = stringr::str_extract(var, '(?<=\\[)[0-9]+(?=\\])')) %>%
    dplyr::mutate(val = basename(gsub('\\\\', '/', val))) %>%
    dplyr::select(file_n = var, sample_file = val) %>%
    dplyr::left_join(tbl_mtd %>%
                       dplyr::filter(grepl("^assay.*ms_run_ref", var)) %>%
                       dplyr::mutate(var = stringr::str_extract(var, '(?<=\\[)[0-9]+(?=\\])')) %>%
                       dplyr::mutate(val = stringr::str_extract(val, '(?<=\\[)[0-9]+(?=\\])')) %>%
                       dplyr::select(assay_id = var, file_n = val),
                     by = 'file_n') %>%
    dplyr::left_join(tbl_mtd %>%
                       dplyr::filter(grepl("^study_var.*assay", var)) %>%
                       tidyr::separate_rows(val, sep="\\,\\s") %>%
                       dplyr::mutate(var = stringr::str_extract(var, '(?<=\\[)[0-9]+(?=\\])')) %>%
                       dplyr::mutate(val = stringr::str_extract(val, '(?<=\\[)[0-9]+(?=\\])')) %>%
                       dplyr::select(studyvar_id = var, assay_id = val),
                     by = 'assay_id') %>%
    dplyr::left_join(tbl_mtd %>%
                       dplyr::filter(grepl("^study_var.*description", var)) %>%
                       dplyr::mutate(var = stringr::str_extract(var, '(?<=\\[)[0-9]+(?=\\])')) %>%
                       dplyr::select(studyvar_id = var, sample_name = val),
                     by = 'studyvar_id') %>%
    as.data.frame()

  # PSM TABLE
  # gather and simplify the PSM data
  tbl_psm <- rfn[which(grepl("^PS[HM]", rfn))] %>%
    stringr::str_replace_all("null", "") %>%
    lapply(vec_tbl) %>%
    dplyr::bind_rows() %>%
    tbl_names(tbl_files) %>%
    dplyr::group_by(sequence, modifications, spectra_ref, charge) %>%
    dplyr::slice_head(n=1) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(sequence, modifications, charge, calc_mass_to_charge) %>%
    dplyr::summarise(num_psms = dplyr::n()) %>%
    dplyr::arrange(dplyr::desc(num_psms)) %>%
    dplyr::rename(mass_to_charge = calc_mass_to_charge)

  # PSM TABLE
  # substitute modification name values
  for( im in 1:nrow(tbl_mods) ){
    id <- tbl_mods$id[im]
    str <- tbl_mods$str[im]
    tbl_psm <- tbl_psm %>% dplyr::mutate(modifications = gsub(id, str, modifications))
  }

  # PSM TABLE
  # deal with modifications missing in PEP that are present in PSM
  # deal with modification site ambiguity
  #
  combine_mods <- function(x){
    if(length(x) <= 1) {return (x)}
    x <- x %>% lapply(stringr::str_split, pattern = ",") %>% lapply(unlist)
    y <- x[[1]]
    z <- c()
    for(i in 2:length(x)){
      z <- c(z, setdiff(y, x[[i]]))
      y <- intersect(x[[i]], y)
      z <- c(z, setdiff(x[[i]], y))
    }
    y <- paste(y, collapse = ",")
    zd <- gsub("[0-9\\-]+", "", z) %>% unique()
    zn <- gsub("[a-zA-Z\\-]+", "", paste(z, collapse = "|"))
    x <- z <- paste0("(", zn, ")-", zd)

    if(y != "") {x <- paste(c(y, z), collapse = ",")}

    return(x)
  }

  tbl_psm <- tbl_psm %>%
    dplyr::group_by(sequence, charge, mass_to_charge) %>%
    dplyr::summarise(
      modifications = combine_mods(modifications),
      num_psms = sum(num_psms),
      .groups = 'drop')

  # PEP TABLE
  # gather PEP data
  # combine PSM data into PEP data
  tbl_pep <- rfn[which(grepl("^PE[HP]", rfn))] %>%
    stringr::str_replace_all("null", "") %>%
    lapply(vec_tbl) %>%
    dplyr::bind_rows() %>%
    tbl_names(tbl_files) %>%
    dplyr::select(!dplyr::matches('modifications|unique')) %>%
    dplyr::left_join(
      tbl_psm, by = c("sequence", "charge", "mass_to_charge")
    )

  # PEP TABLE
  # combine PSM data into PEP data
  tbl_pep <- tbl_pep %>%
    dplyr::group_by(sequence, accession, modifications) %>%
    dplyr::summarise(num_psms = sum(num_psms),
                     .groups = 'drop') %>%
    dplyr::left_join(
      tbl_pep %>%
        dplyr::select(!dplyr::matches('num\\_')) %>%
        dplyr::group_by(sequence, accession, modifications) %>%
        dplyr::slice_head(n = 1) %>%
        dplyr::ungroup(),
      by = c('sequence', 'accession', 'modifications')
    )

  # PRT TABLE
  # gather PRT data
  tbl_pro <- rfn[which(grepl("^PR[HT]", rfn))] %>%
    stringr::str_replace_all("null", "") %>%
    lapply(vec_tbl) %>%
    dplyr::bind_rows() %>%
    tbl_names(tbl_files)

  if(analyte == 'peptides') {
    tbl_dat <- tbl_pep %>%
      dplyr::left_join(tbl_pro %>%
                         dplyr::select(dplyr::matches("^accession|description$")),
                       by = 'accession')
  } else {
    # collapse the merged tables down to protein-level data
    tbl_dat <- tbl_pro %>%
      dplyr::left_join(tbl_pep %>%
                         dplyr::group_by(accession) %>%
                         dplyr::summarise(
                           num_peptides = dplyr::n(),
                           num_unique_peptides = length(unique(sequence)),
                           num_psms = sum(num_psms),
                           .groups = 'drop'
                         ), by = "accession")
  }

  return(list(platform = platform, analyte = analyte, data = tbl_dat))
}
