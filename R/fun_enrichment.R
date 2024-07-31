#' Helper function to calculate term enrichment
#'
#' @param data tidyproteomics data table object
#' @param x the annotation to compute enrichment for
#'
#' @return list of vectors
#'
calc_enrichment <- function(
    data,
    x
){
  stats::median(data$log2_foldchange[which(data$annotation == x)], na.rm = T) /
    stats::median(data$log2_foldchange, na.rm = T)
}
