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
  abundance_raw <- NULL
  imputed <- NULL

  if(is.null(identifier)) {cli::cli_abort("No identifier provided for codifying")}
  if(is.null(table)) {cli::cli_abort("No table provided for codifying")}
  if(length(intersect(annotations, names(table))) == 0) { annotations <- NULL }

  cols_experiments <- c("import_file","sample_file","sample_id","sample","replicate")

  get_experiments <- paste(cols_experiments, collapse="|")
  get_quantitative <- paste(c(paste0("^", identifier, "$"), "^sample$", "replicate", "sample_id", "abundance"), collapse="|")
  get_annotations <- paste0("^", paste(unique(c(identifier, annotations)), collapse="$|^"), "$")
  get_accounting <- paste(c(identifier, "sample_id", "impute", "match_between", "num\\_", "\\_group"), collapse="|")

  # remove the replicate column to avoid merge errors
  table[,'replicate'] <- NULL

  tb_experiments <- table %>%
    dplyr::select(dplyr::matches(get_experiments)) %>%
    unique() %>%
    dplyr::group_by(sample) %>%
    dplyr::arrange(sample_id) %>%
    dplyr::mutate(replicate = as.character(dplyr::row_number())) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(sample, replicate)

  table <- table %>%
    tibble::as_tibble() %>%
    dplyr::full_join(tb_experiments,
                     by = intersect(names(table), cols_experiments))

  tb_experiments <- table %>%
    dplyr::select(dplyr::matches(get_experiments)) %>%
    unique() %>%
    dplyr::relocate('sample_id')

  tb_quantitative <- table %>%
    dplyr::select(dplyr::matches(get_quantitative)) %>%
    dplyr::select(!dplyr::matches("num\\_|\\_group")) %>%
    dplyr::relocate(c('sample_id','sample','replicate'))

  tb_quantitative_tmp <- table %>%
    dplyr::group_by_at(setdiff(colnames(tb_quantitative), 'abundance_raw')) %>%
    dplyr::summarise(
      abundance_raw = sum(abundance_raw),
      .groups = 'drop'
    )

  if(nrow(tb_quantitative) > nrow(tb_quantitative_tmp)){
    cli::cli_alert_info(c("i" = "... peptide accounting indicates multiple precursors"))
    cli::cli_alert_info(c("i" = "...... {.emph {nrow(tb_quantitative)}} quantitative values merged down to {.emph {nrow(tb_quantitative_tmp)}}"))
    tb_quantitative <- tb_quantitative_tmp
  }

  if(!is.null(annotations) & length(annotations) > 0) {
    tb_annotations <- table %>%
      dplyr::select(dplyr::matches(get_annotations))

    if(ncol(tb_annotations) == 1) {
      tb_annotations <- c()
    } else {
      tb_annotations <- tb_annotations %>%
        dplyr::select(!dplyr::matches("num\\_|\\_group")) %>%
        tidyr::pivot_longer(
          cols = !dplyr::matches(paste(paste0("^", identifier, "$"), collapse = "|")),
          names_to = 'term',
          values_to = 'annotation'
        ) %>% unique() %>%
        dplyr::filter(!is.na(annotation))
    }
  } else {
    tb_annotations <- c()
  }

  tb_accounting <- table %>%
    dplyr::select(dplyr::matches(get_accounting)) %>%
    dplyr::relocate("sample_id") %>%
    unique()

  # need to understand what columns are in accounting
  cols_share <- intersect(colnames(tb_accounting), colnames(tb_quantitative))
  cols_grouped <- setdiff(colnames(tb_accounting), cols_share)

  # test and indicate that imputation was not accounted for
  if(length(intersect(c("match_between_runs","imputed"), names(tb_accounting))) == 0) {
    tb_accounting <- tb_accounting %>% dplyr::mutate(imputed = 0)
  }

  # merge MBR and imputed into a single is.imputed value
  tb_accounting <- tb_accounting %>%
    tidyr::pivot_longer(
      dplyr::matches("match_between_runs|imputed"),
      names_to = 'method',
      values_to = 'imputed'
    )  %>%
    dplyr::mutate(imputed = as.numeric(imputed)) %>%
    dplyr::select(!c('method')) %>%
    dplyr::group_by_at(cols_share) %>%
    dplyr::slice_max(imputed, n=1, with_ties = FALSE) %>%
    dplyr::ungroup()

  n_rows_tb_accounting <- tb_accounting %>% dplyr::select(c('sample_id', identifier)) %>% unique() %>% nrow()
  n_rows_tb_quantitative <- tb_quantitative %>% dplyr::select(c('sample_id', identifier)) %>% unique() %>% nrow()

  if(n_rows_tb_accounting != n_rows_tb_quantitative) {
    cli::cli_div(theme = list(span.emph = list(color = "#ff4500")))
    cli::cli_alert_info(c("x" = "Something did not parse correctly"))
    cli::cli_alert_info(c("i" = "... quantitative values shows {.emph {nrow(tb_quantitative)}} unique values"))
    cli::cli_alert_info(c("i" = "... accounting values shows {.emph {nrow(tb_accounting)}} unique values"))
    cli::cli_alert_info(c("i" = "... check the sample group expressions for {.emph pattern_extract}"))
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
