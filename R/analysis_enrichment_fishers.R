#' A function for evaluating term enrichment via Fischer's Exact method
#'
#' @param data_expression a tibble from and two sample expression difference analysis
#' @param data tidyproteomics data object
#' @param term_group a character string referencing "term" in the annotations table
#' @param log2fc_min a numeric defining the minimum log2 foldchange to consider as "enriched"
#' @param significance_max a numeric defining the maximum statistical significance to consider as "enriched"
#' @param cpu_cores the number of threads used to speed the calculation
#' @param ... pass through arguments
#'
#' @return a tibble
#'
enrichment_fishersexact <- function(
    data_expression = NULL,
    data = NULL,
    term_group = NULL,
    log2fc_min = 0,
    significance_min = 0.05,
    cpu_cores = 1,
    ...
){

  # visible bindings
  log2_foldchange <- NULL
  annotation <- NULL
  matches <- NULL
  p.value <- NULL
  p_value <- NULL
  adj_p_value <- NULL

  term_group <- rlang::arg_match(term_group, get_annotation_terms(data))
  check_data(data)

  identifier <- data$identifier
  tbl_x <- data_expression %>%
    dplyr::arrange(dplyr::desc(log2_foldchange)) %>%
    dplyr::mutate(rank = dplyr::row_number()) %>%
    dplyr::left_join(get_annotations(data, term_group),
                     by = data$identifier) %>%
    dplyr::mutate(annotation = ifelse(is.na(annotation), 'other', annotation)) %>%
    dplyr::mutate(annotation = as.character(annotation))

  colnames(tbl_x)[which(colnames(tbl_x) == data$identifier)] <- 'identifier'

  # calculating p-values from gene set enrichment
  # for over-representation analysis
  # ----------------------------------------------
  #                  :    in    :  not in  :
  #                  : gene set : gene set : TOTAL
  # ----------------------------------------------
  #  differential    :          :          :
  #  expressed       :    d     :    b     :   k
  # ----------------------------------------------
  #  NOT differential:          :          :
  #  expressed       :    q     :    p     :   -
  # ----------------------------------------------
  #  TOTAL           :    m     :    n     :   t
  # ----------------------------------------------

  list.append <- function (x, i){x[[length(x) + 1]] <- i; x}
  enrichment <- function(data, x){
    stats::median(data$log2_foldchange[which(data$annotation == x)], na.rm = T) /
      stats::median(data$log2_foldchange, na.rm = T)
  }

  tbl_x_sig <- tbl_x %>%
    dplyr::filter(p_value <= significance_min) %>%
    dplyr::filter(log2_foldchange >= log2fc_min)

  tbl_out <- list()
  annotations <- tbl_x %>% dplyr::filter(annotation != 'other') %>% dplyr::select(annotation) %>% unlist() %>% unique()
  for(test_annotation in annotations){

    # cli::cli_alert_info("{test_annotation}")

    tbl_x_term <- tbl_x %>% dplyr::filter(annotation == test_annotation)
    tbl_x_term_enrichment <- tbl_x %>% calc_enrichment(test_annotation)
    tbl_x_sig_term <- tbl_x_sig %>% dplyr::filter(annotation == test_annotation)

    d <- tbl_x_sig_term$identifier %>% unique() %>% length()
    m <- tbl_x_term$identifier %>% unique() %>% length()
    q <- m - d
    k <- tbl_x_sig$identifier %>% unique() %>% length()
    t <- tbl_x$identifier %>% unique() %>% length()
    b <- k - d
    n <- t - m
    p <- n - b

    # hypergeometric
    # pval_h <- 1 - phyper(d - 1, m, n, k)
    # cli::cli_alert_info("\t hypergeometric {pval}")

    # binomial
    # pval_b <- 1 - pbinom(d - 1, m, k / t)
    # cli::cli_alert_info("\t binomial {pval}")

    # Fisher's Exact
    cm = matrix(c(  d, q, b, p  ), nrow = 2)
    pval <- fisher.test(cm)
    # cli::cli_alert_info("\t Fisher's Exact {pval$p.value}")

    tbl_out <- tbl_out %>% list.append(tibble::tibble(annotation = test_annotation,
                                                     p_value = pval$p.value,
                                                     adj_p_value = NA,
                                                     enrichment = tbl_x_term_enrichment,
                                                     size = m))
  }

  tbl_out <- tbl_out %>%
    dplyr::bind_rows() %>%
    dplyr::arrange(p_value) %>%
    dplyr::mutate(adj_p_value = p.adjust(p_value))

  return(tbl_out)

}
