#' A function for evaluating term enrichment via Wilcoxon Rank Sum
#'
#' @param data_expression a tibble from and two sample expression difference analysis
#' @param data tidyproteomics data object
#' @param term_group a character string referencing "term" in the annotations table
#' @param cpu_cores the number of threads used to speed the calculation
#' @param ... pass through arguments
#'
#' @return a tibble
#'
enrichment_wilcoxon <- function(
    data_expression = NULL,
    data = NULL,
    term_group = NULL,
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
  tbl_expression <- data_expression %>%
    dplyr::arrange(dplyr::desc(log2_foldchange)) %>%
    dplyr::mutate(rank = dplyr::row_number()) %>%
    dplyr::left_join(get_annotations(data, term_group),
                     by = data$identifier) %>%
    dplyr::mutate(annotation = ifelse(is.na(annotation), 'other', annotation)) %>%
    dplyr::mutate(annotation = as.character(annotation))

  colnames(tbl_expression)[which(colnames(tbl_expression) == data$identifier)] <- 'identifier'

  # calculate significance by wilcoxon rank
  wilcox_test <- function(data, x){
    stats::wilcox.test(
      data$rank[which(data$annotation == x)],
      data$rank[which(data$annotation != x)],
    )
  }

  f_enrich <- function(annotation_str, x){

    tryCatch({

      # table of annotation to test
      tbl_test <- x %>%
        dplyr::group_by(identifier) %>%
        dplyr::mutate(rank = ifelse(annotation == annotation_str, 0, rank)) %>%
        dplyr::mutate(annotation = ifelse(annotation != annotation_str, 'other', annotation)) %>%
        dplyr::slice_min(rank, n = 1, with_ties = FALSE) %>%
        dplyr::ungroup() %>%
        dplyr::arrange(dplyr::desc(log2_foldchange)) %>%
        dplyr::mutate(rank = dplyr::row_number())

      out <- tbl_test %>%
        wilcox_test(annotation_str) %>%
        broom::tidy() %>%
        dplyr::mutate(enrichment = tbl_test %>% calc_enrichment(annotation_str),
                      annotation = annotation_str,
                      size = tbl_test %>%
                        dplyr::filter(annotation == annotation_str) %>%
                        nrow()) %>%
        dplyr::rename(p_value = p.value) %>%
        dplyr::select(!matches('statistic|method|alternative')) %>%
        dplyr::relocate(annotation)

      return(out)

    }, error = function(err) {
      cli::cli_div(theme = list(span.emph = list(color = "#ff4500")))
      cli::cli_alert_info("annotation {.emph {annotation_str}} had issues, not reported")

      return(NULL)
    })

  }

  data_out <- unique(tbl_expression$annotation) %>%
    parallel::mclapply(f_enrich, tbl_expression, mc.cores = cpu_cores) %>%
    dplyr::bind_rows() %>%
    dplyr::mutate(adj_p_value = stats::p.adjust(p_value)) %>%
    dplyr::relocate(adj_p_value, .after = p_value) %>%
    dplyr::arrange(p_value) %>%
    dplyr::filter(size >= 3) %>%
    dplyr::filter(size <= length(unique(tbl_expression$identifier))*.66)

  return(data_out)
}
