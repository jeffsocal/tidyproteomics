#' @importFrom prettyunits pretty_bytes
#' @importFrom stats fisher.test p.adjust var
#' @importFrom utils head
NULL

# Global variable bindings to silence R CMD check codetools notes

utils::globalVariables(c(
  ".data", "DF", "Intensity", "MissingPercentage", "Protein", "SE", "Tvalue",
  "abundance", "adj.pvalue", "adj_p_value", "annotation", "average_expression",
  "enrichment", "foldchange", "identifier", "imputed", "is_enriched",
  "is_sig", "is_significant", "keep", "label_pos", "log2FC", "log2_abundance",
  "log2_foldchange", "match_between_runs", "matches", "max_y", "metric",
  "min_group", "min_y", "n", "origin", "p_value", "pvalue",
  "proportional_expression", "replicate", "sample", "sample_file", "sample_id",
  "sample_origin", "size", "term", "var", "y_bar", "y_span", "y_text"
))
