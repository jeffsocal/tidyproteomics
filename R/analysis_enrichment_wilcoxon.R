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
  p_value <- NULL
  adj_p_value <- NULL

  term_group <- rlang::arg_match(term_group, get_annotation_terms(data))
  check_data(data)

  id_col <- data$identifier[1]

  # 1. Clean expression data: remove NA / Inf and deduplicate identifiers
  df_stats <- data_expression %>%
    dplyr::select(identifier = dplyr::all_of(id_col), log2_foldchange) %>%
    dplyr::filter(!is.na(log2_foldchange), !is.infinite(log2_foldchange)) %>%
    dplyr::distinct(identifier, .keep_all = TRUE) %>%
    dplyr::arrange(dplyr::desc(log2_foldchange))

  if(nrow(df_stats) == 0) {
    return(tibble::tibble(
      annotation = character(0),
      p_value = numeric(0),
      adj_p_value = numeric(0),
      enrichment = numeric(0),
      size = integer(0)
    ))
  }

  # 2. Compute ranks once for the whole experiment (average ranks handle ties cleanly)
  df_stats$rank <- rank(-df_stats$log2_foldchange, ties.method = "average")
  c_stats <- stats::setNames(df_stats$log2_foldchange, df_stats$identifier)
  c_ranks <- stats::setNames(df_stats$rank, df_stats$identifier)
  all_ids <- names(c_ranks)
  n_total <- length(c_ranks)

  # 3. Clean annotations: remove NA, blank, and 'other'
  tbl_anno <- get_annotations(data, term_group) %>%
    dplyr::filter(!is.na(annotation), annotation != "", annotation != "other")

  if(nrow(tbl_anno) == 0) {
    return(tibble::tibble(
      annotation = character(0),
      p_value = numeric(0),
      adj_p_value = numeric(0),
      enrichment = numeric(0),
      size = integer(0)
    ))
  }

  # 4. Map pathways and pre-filter by size BEFORE testing
  pathways <- split(as.character(tbl_anno[[id_col]]), tbl_anno$annotation)
  pathways <- lapply(pathways, function(p) intersect(unique(p), all_ids))

  sizes <- lengths(pathways)
  keep <- sizes >= 3 & sizes <= n_total * 0.66
  pathways <- pathways[keep]
  sizes <- sizes[keep]

  if(length(pathways) == 0) {
    return(tibble::tibble(
      annotation = character(0),
      p_value = numeric(0),
      adj_p_value = numeric(0),
      enrichment = numeric(0),
      size = integer(0)
    ))
  }

  # 5. Fast Wilcoxon test over clean vectors
  f_test <- function(term_name) {
    p <- pathways[[term_name]]
    in_ranks <- c_ranks[p]
    out_ranks <- c_ranks[!(all_ids %in% p)]

    tryCatch({
      wt <- stats::wilcox.test(in_ranks, out_ranks, exact = FALSE)
      wt$p.value
    }, error = function(e) {
      NA_real_
    })
  }

  pvals <- if (cpu_cores > 1 && .Platform$OS.type != "windows") {
    res_list <- parallel::mclapply(names(pathways), f_test, mc.cores = cpu_cores)
    vapply(res_list, function(v) if(is.numeric(v) && length(v) == 1) v else NA_real_, numeric(1))
  } else {
    vapply(names(pathways), f_test, numeric(1))
  }

  med_enrichment <- vapply(pathways, function(p) {
    stats::median(c_stats[p], na.rm = TRUE)
  }, numeric(1))

  data_out <- tibble::tibble(
    annotation = names(pathways),
    p_value = unname(pvals),
    adj_p_value = stats::p.adjust(unname(pvals), method = "BH"),
    enrichment = unname(med_enrichment),
    size = unname(sizes)
  ) %>%
    dplyr::filter(!is.na(p_value)) %>%
    dplyr::arrange(p_value)

  return(data_out)
}
