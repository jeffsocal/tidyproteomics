#' A function for evaluating term enrichment via GSEA
#'
#' @param data_expression a tibble from and two sample expression difference analysis
#' @param data tidyproteomics data object
#' @param term_group a character string referencing "term" in the annotations table
#'
#' @return a tibble
#'
#'
enrichment_gsea <- function(
    data_expression = NULL,
    data = NULL,
    term_group = NULL
){

  # visible bindings
  desc <- NULL
  log2_foldchange <- NULL
  annotation <- NULL
  pval <- NULL
  pathway <- NULL
  padj <- NULL
  ES <- NULL
  NES <- NULL
  log2err <- NULL
  p_value <- NULL

  term_group <- rlang::arg_match(term_group, get_annotation_terms(data))
  check_data(data)

  data_expression <- data_expression %>%
    munge_identifier(identifiers = data$identifier) %>%
    dplyr::inner_join(get_annotations(data, term_group) %>%
                        munge_identifier(identifiers = data$identifier),
                      by = 'identifier') %>%
    dplyr::arrange(desc(log2_foldchange)) %>%
    dplyr::mutate(rank = dplyr::row_number(),
                  nrank = (1 - rank / max(rank)) * 2)

  data_term <- data_expression %>%
    tidyr::separate_rows(annotation, sep=";") %>%
    dplyr::select(c('identifier', 'annotation')) %>%
    dplyr::group_by(annotation) %>%
    tidyr::nest()

  l_term <- list()
  for( i in 1:nrow(data_term) ){
    if(is.na(data_term$annotation[i] %>% unlist()))
      next
    l_term[[data_term$annotation[i] %>% unlist()]] <- data_term$data[i] %>% unlist() %>% as.character()
  }

  c_rank <- data_expression$rank %>% unlist() %>% as.numeric()
  names(c_rank) <- data_expression$identifier %>% unlist() %>% as.character()

  tryCatch({
    data_out <- fgsea::fgsea(pathways = l_term,
                             stats = c_rank,
                             scoreType = 'pos',
                             minSize=15,
                             maxSize=600) %>%
      tibble::as_tibble() %>%
      dplyr::arrange(pval) %>%
      dplyr::select(
        annotation = pathway,
        p_value = pval,
        adj_p_value = padj,
        enrichment = ES,
        enrichment_normalized = NES,
        log2err
      ) %>%
      dplyr::full_join(
        data_expression %>%
          dplyr::group_by(annotation) %>%
          dplyr::summarise(size = dplyr::n(),
                           .groups = 'drop'),
        by = c('annotation')
      )

  }, error = function(err) {
    cli::cli_alert_info("{annotation_str}: had issues")
  })

  data_out <- data_out %>%
    dplyr::bind_rows() %>%
    dplyr::arrange(p_value)

  return(data_out)
}
