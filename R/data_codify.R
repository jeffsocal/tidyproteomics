#' Build a tidyproteomics data object
#'
#' @description
#' `data_codify()` is a helper function
#'
#' @param table tidyproteomics data object
#' @param identifier a character vector
#' @param annotations a character vector
#'
#' @return tidyproteomics data object
#'
codify <- function(
    table=NULL,
    identifier=NULL,
    annotations = NULL
){

  # visible bindings
  import_file <- NULL
  sample_file <- NULL
  annotation <- NULL

  if(is.null(identifier)) {cli::cli_abort("No identifier provided for codifying")}
  if(is.null(table)) {cli::cli_abort("No table provided for codifying")}
  if(length(intersect(annotations, names(table))) == 0) { annotations <- NULL }

  get_experiments <- paste(c("import_file","sample_file","sample_id","sample","replicate"), sep="|")
  get_quantitative <- paste(c(paste0("^", identifier, "$"), "^sample$", "replicate", "sample_id", "abundance"), sep="|")
  get_annotations <- paste(c(identifier, annotations), sep="|")
  get_accounting <- paste(c(identifier, "sample_id", "impute", "match_between", "num\\_", "\\_group"), sep="|")

  tb_experiments = table %>%
    dplyr::select(dplyr::matches(get_experiments)) %>%
    unique() %>%
    dplyr::arrange(import_file, sample_file) %>%
    dplyr::group_by(sample) %>%
    dplyr::mutate(replicate = as.character(dplyr::row_number())) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(sample, replicate)

  table <- table %>%
    tibble::as_tibble() %>%
    dplyr::full_join(tb_experiments,
                     by = intersect(names(table), get_experiments))

  tb_experiments <- table %>%
    dplyr::select(dplyr::matches(get_experiments)) %>%
    unique() %>%
    dplyr::relocate('sample_id')

  tb_quantitative <- table %>%
    dplyr::select(dplyr::matches(get_quantitative)) %>%
    dplyr::select(!dplyr::matches("num\\_|\\_group")) %>%
    dplyr::relocate(c('sample_id','sample','replicate'))

  if(!is.null(annotations) || length(annotations) > 0) {
    tb_annotations <- table %>% dplyr::select(dplyr::matches(get_annotations)) %>%
      dplyr::select(!dplyr::matches("num\\_|\\_group")) %>%
      tidyr::pivot_longer(
        cols = !dplyr::matches(paste(paste0("^", identifier, "$"), collapse = "|")),
        names_to = 'term',
        values_to = 'annotation'
      ) %>% unique() %>%
      dplyr::filter(!is.na(annotation))
  } else {
    tb_annotations <- c()
  }

  tb_accounting <- table %>%
    dplyr::select(dplyr::matches(get_accounting)) %>%
    dplyr::select(!dplyr::matches("\\_name$")) %>%
    dplyr::relocate("sample_id") %>%
    unique()

  if(nrow(tb_accounting) != nrow(tb_quantitative)) {
    cli::cli_div(theme = list(span.emph = list(color = "#ff4500")))
    cli::cli_alert_info(c("x" = "Something did not parse correctly",
                     "i" = "Table of quantitative values has {.emph {nrow(tb_quantitative)}} rows",
                     "i" = "Table of accounting values has {.emph {nrow(tb_accounting)}} rows",
                     "i" = "Check the sample group expressions for {.emph pattern_extract}"))
  }

  list(
    experiments = tb_experiments,
    quantitative = tb_quantitative,
    accounting = tb_accounting,
    annotations = tb_annotations
  )
}

#' Meld a tidyproteomics data object into a single table
#'
#' @description
#' `data_meld()` is a helper function
#'
#' @param data tidyproteomics data object
#' @param single_quant_source a boolean to indicate if only a single quantitative value should be reported
#'
#' @return a tibble
#'
meld <- function(
    data=NULL,
    single_quant_source = FALSE
){

  # visible bindings
  annotation <- NULL
  term <- NULL

  check_data(data)

  identifier <- data$identifier
  join_second <- c(identifier, "sample_id", "sample", "replicate")

  tb_quant <- data$quantitative

  if(single_quant_source == TRUE) {
    # select only a single quantitative value to move forward
    quant_source <- data$quantitative_source
    tb_quant <- tb_quant %>%
      dplyr::select(!tidyselect::matches('abundance') | tidyselect::matches(paste0(quant_source, "$"))) %>%
      dplyr::rename(abundance = tidyselect::matches(quant_source))
  }

  tb_accnt <- data$accounting %>%
    dplyr::left_join(data$experiments, by = 'sample_id')

  if(is.null(data$annotations) || nrow(data$annotations) == 0) {
  } else {
    tb_annot <- data$annotations %>%
      dplyr::filter(!is.na(term), !is.na(annotation)) %>%
      tidyr::pivot_wider(names_from = 'term', values_from = 'annotation')
    join_first <- intersect(identifier, names(tb_annot))

    tb_quant <- tb_quant %>%
      dplyr::left_join(tb_annot, by = join_first)
  }

  tb_quant <- tb_quant %>%
    dplyr::left_join(tb_accnt, by = join_second)


  return(tb_quant)
}
