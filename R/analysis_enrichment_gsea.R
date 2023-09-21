#' A function for evaluating term enrichment via GSEA
#'
#' @param data_expression a tibble from and two sample expression difference analysis
#' @param data tidyproteomics data object
#' @param term_group a character string referencing "term" in the annotations table
#' @param score_type a character string used in the fgsea package
#' @param cpu_cores the number of threads used to speed the calculation
#'
#' @return a tibble
#'
#'
enrichment_gsea <- function(
    data_expression = NULL,
    data = NULL,
    term_group = NULL,
    score_type = c("std", "pos", "neg"),
    cpu_cores = 1
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

  tbl_expression <- data_expression %>%
    dplyr::arrange(dplyr::desc(log2_foldchange)) %>%
    dplyr::mutate(rank = dplyr::row_number()) %>%
    dplyr::left_join(get_annotations(data, term_group),
                     by = data$identifier) %>%
    dplyr::mutate(annotation = ifelse(is.na(annotation), 'other', annotation)) %>%
    dplyr::mutate(annotation = as.character(annotation))

  colnames(tbl_expression)[which(colnames(tbl_expression) == data$identifier)] <- 'identifier'

  f_enrich <- function(annotation_str, x){

    # table of annotation to test
    tbl_test <- x %>%
      dplyr::group_by(identifier) %>%
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
      out <- fgsea::fgsea(pathways = l_term,
                          stats = c_rank,
                          scoreType = score_type,
                          minSize=3,
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
        dplyr::filter(annotation != 'other') %>%
        dplyr::inner_join(
          tbl_test %>%
            dplyr::group_by(annotation) %>%
            dplyr::summarise(size = dplyr::n(),
                             .groups = 'drop'),
          by = c('annotation')
        )

      return(out)

    }, error = function(err) {
      cli::cli_div(theme = list(span.emph = list(color = "#ff4500")))
      cli::cli_alert_info("annotation {.emph {annotation_str}} had issues, not reported")
      cli::cli_abort(err)

      return(NULL)
    })

  }

  data_out <- unique(tbl_expression$annotation) %>%
    parallel::mclapply(f_enrich, tbl_expression, mc.cores = cpu_cores) %>%
    dplyr::bind_rows() %>%
    dplyr::mutate(adj_p_value = stats::p.adjust(p_value)) %>%
    dplyr::relocate(adj_p_value, .after = p_value) %>%
    dplyr::arrange(p_value)

  return(data_out)
}
