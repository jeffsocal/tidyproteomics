#' A function for evaluating term enrichment via GSEA
#'
#' @param data_expression a tibble from and two sample expression difference analysis
#' @param data tidyproteomics data object
#' @param term_group a character string referencing "term" in the annotations table
#' @param score_type a character string.
#'
#' @return a tibble
#'
#'
enrichment_gsea <- function(
    data_expression = NULL,
    data = NULL,
    term_group = NULL,
    score_type = c("std", "pos", "neg")
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
  identifier <- NULL
  adj_p_value <- NULL

  set.seed(1234)

  term_group <- rlang::arg_match(term_group, get_annotation_terms(data))
  score_type <- rlang::arg_match(score_type)
  check_data(data)

  data_expression <- data_expression %>%
    dplyr::left_join(get_annotations(data, term_group),
                     by = data$identifier)

  data_out <- list()
  for(annotation_str in unique(data_expression$annotation)) {

    if(is.na(annotation_str)) { next }

      # table of annotation to test
      tbl_test <- data_expression %>%
        tidyr::unite(identifier, data$identifier) %>%
        dplyr::group_by(identifier) %>%
        dplyr::mutate(rank = dplyr::row_number()) %>%
        dplyr::mutate(annotation = paste0("", as.character(annotation))) %>%
        dplyr::mutate(rank = ifelse(annotation == annotation_str, 0, rank)) %>%
        dplyr::mutate(annotation = ifelse(annotation != annotation_str, 'other', annotation)) %>%
        dplyr::slice_min(rank, n = 1, with_ties = FALSE) %>%
        dplyr::ungroup() %>%
        dplyr::arrange(dplyr::desc(log2_foldchange)) %>%
        dplyr::mutate(rank = dplyr::row_number())

    tryCatch({

      c_rank <- tbl_test$rank %>% unlist() %>% as.numeric()
      names(c_rank) <- tbl_test$identifier %>% unlist() %>% as.character()

      data_term <- tbl_test %>%
        dplyr::select(c('identifier', 'annotation')) %>%
        dplyr::group_by(annotation) %>%
        tidyr::nest()

      l_term <- list()
      for( i in 1:nrow(data_term) ){
        if(is.na(data_term$annotation[i] %>% unlist()))
          next
        l_term[[data_term$annotation[i] %>% unlist()]] <- data_term$data[i] %>% unlist() %>% as.character()
      }

      # fgsea is set up to run all possible terms, assuming only one identifier per term
      data_out[[annotation_str]] <- fgsea::fgsea(pathways = l_term,
                               stats = c_rank,
                               scoreType = score_type,
                               minSize=15,
                               maxSize=nrow(tbl_test)*.66) %>%
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
          tbl_test %>%
            dplyr::group_by(annotation) %>%
            dplyr::summarise(size = dplyr::n(),
                             .groups = 'drop'),
          by = c('annotation')
        )


    }, error = function(err) {
      cli::cli_div(theme = list(span.emph = list(color = "#ff4500")))
      cli::cli_alert_info("annotation {.emph {annotation_str}} had issues, not reported")
      cli::cli_abort(err)
    })
  }

  data_out <- data_out %>%
    dplyr::bind_rows() %>%
    dplyr::filter(annotation != 'other') %>%
    dplyr::mutate(adj_p_value = stats::p.adjust(p_value)) %>%
    dplyr::relocate(adj_p_value, .after = p_value) %>%
    dplyr::arrange(p_value)

  return(data_out)
}
