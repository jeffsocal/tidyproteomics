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
  annotation <- NULL
  pval <- NULL
  pathway <- NULL
  padj <- NULL
  ES <- NULL
  NES <- NULL
  log2err <- NULL
  p_value <- NULL
  adj_p_value <- NULL
  size <- NULL

  term_group <- rlang::arg_match(term_group, get_annotation_terms(data))
  score_type <- rlang::arg_match(score_type)
  check_data(data)

  id_col <- data$identifier[1]
  c_stats <- data_expression$log2_foldchange
  names(c_stats) <- as.character(data_expression[[id_col]])
  c_stats <- c_stats[!is.na(c_stats) & !is.infinite(c_stats)]
  c_stats <- sort(c_stats, decreasing = TRUE)

  tbl_anno <- get_annotations(data, term_group)
  tbl_anno <- tbl_anno %>% dplyr::filter(!is.na(annotation), annotation != "", annotation != "other")

  pathways <- split(as.character(tbl_anno[[id_col]]), tbl_anno$annotation)
  pathways <- lapply(pathways, function(p) intersect(unique(p), names(c_stats)))

  bpparam <- if (cpu_cores > 1 && .Platform$OS.type != "windows") {
    BiocParallel::MulticoreParam(workers = cpu_cores)
  } else {
    BiocParallel::SerialParam()
  }

  out <- tryCatch({
    fgsea::fgsea(
      pathways = pathways,
      stats = c_stats,
      scoreType = score_type,
      minSize = 3,
      maxSize = length(c_stats) * 0.66,
      BPPARAM = bpparam
    )
  }, error = function(err) {
    fgsea::fgsea(
      pathways = pathways,
      stats = c_stats,
      scoreType = score_type,
      minSize = 3,
      maxSize = length(c_stats) * 0.66,
      BPPARAM = BiocParallel::SerialParam()
    )
  })

  data_out <- out %>%
    tibble::as_tibble() %>%
    dplyr::select(
      annotation = pathway,
      p_value = pval,
      adj_p_value = padj,
      enrichment = ES,
      enrichment_normalized = NES,
      log2err,
      size
    ) %>%
    dplyr::arrange(p_value)

  return(data_out)
}
