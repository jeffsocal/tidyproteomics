#' A function for evaluating term enrichment via Wilcoxon Rank Sum
#'
#' @param data_expression a tibble from and two sample expression difference analysis
#' @param data tidyproteomics data object
#' @param term_group a character string referencing "term" in the annotations table
#' @param ... pass through arguments
#'
#' @return a tibble
#'
enrichment_wilcoxon <- function(
    data_expression = NULL,
    data = NULL,
    term_group = NULL,
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
  data_expression <- data_expression %>%
    dplyr::left_join(get_annotations(data, term_group),
                     by = data$identifier)

  # score = ( (s - p) * sqrt(m) ) / q
  # calculate significance by wilcoxon rank
  wilcox_test <- function(data, x){
    stats::wilcox.test(
      data$rank[which(data$annotation == x)],
      data$rank[which(data$annotation != x)],
    )
  }

  enrichment <- function(data, x){
    stats::median(data$log2_foldchange[which(data$annotation == x)], na.rm = T) /
      stats::median(data$log2_foldchange, na.rm = T)
  }

  data_out <- list()
  for(annotation_str in unique(data_expression$annotation)) {

    if(is.na(annotation_str)) { next }

    # table of annotation to test
    tbl_test <- data_expression %>%
      dplyr::group_by(dplyr::across(identifier)) %>%
      dplyr::mutate(annotation = as.character(annotation)) %>%
      dplyr::mutate(rank = dplyr::row_number()) %>%
      dplyr::mutate(rank = ifelse(annotation == annotation_str, 0, rank)) %>%
      dplyr::slice_min(rank, n = 1, with_ties = FALSE) %>%
      dplyr::ungroup() %>%
      dplyr::arrange(dplyr::desc(log2_foldchange)) %>%
      dplyr::mutate(rank = dplyr::row_number())

    tryCatch({

      data_out[[annotation_str]] <- tbl_test %>%
        wilcox_test(annotation_str) %>%
        broom::tidy() %>%
        dplyr::mutate(enrichment = tbl_test %>% enrichment(annotation_str),
                      annotation = annotation_str,
                      size = data_expression %>%
                        dplyr::filter(annotation == annotation_str) %>%
                        nrow()) %>%
        dplyr::rename(p_value = p.value) %>%
        dplyr::select(!matches('statistic|method|alternative')) %>%
        dplyr::relocate(annotation)
    }, error = function(err) {
      cli::cli_div(theme = list(span.emph = list(color = "#ff4500")))
      cli::cli_alert_info("annotation {.emph {annotation_str}} had issues, not reported")
    })
  }

  data_out <- data_out %>%
    dplyr::bind_rows() %>%
    dplyr::mutate(adj_p_value = stats::p.adjust(p_value)) %>%
    dplyr::arrange(p_value) %>%
    dplyr::relocate(adj_p_value, .after = p_value)

  return(data_out)
}
