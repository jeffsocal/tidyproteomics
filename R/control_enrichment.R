#' Compute protein enrichment
#'
#' @description
#' `enrichment()` is an analysis function that computes the protein summary
#' statistics for a given tidyproteomics data object.
#'
#' @param data tidyproteomics data object
#' @param ... two sample comparison e.g. experimental/control
#' @param .term a character string referencing ".term" in the annotations table
#' @param .method a character string
#'
#' @return a tibble
#' @export
#'
#' @examples
#' library(dplyr, warn.conflicts = FALSE)
#' library(tidyproteomics)
#'
#' # using the default GSEA method
#' ecoli_proteins %>%
#'    expression(ko/wt) %>%
#'    enrichment(ko/wt, .term = 'biological_process') %>%
#'    export_analysis(ko/wt, .analysis = "enrichment")
#'
#' # using a Wilcoxon Rank Sum method
#' ecoli_proteins %>%
#'    expression(ko/wt) %>%
#'    enrichment(ko/wt, .term = 'biological_process', .method = 'wilcoxon') %>%
#'    export_analysis(ko/wt, .analysis = "enrichment")
#'
enrichment <- function(
    data = NULL,
    ...,
    .term = NULL,
    .method = c('gsea', 'wilcoxon')
){

  check_data(data)
  str_quo <- tidyproteomics_quo(...)
  if(is.null(str_quo)) { return(data) }
  experiment <- str_quo['variable']
  control <- str_quo['value']
  if(str_quo['operator'] != "/") {
    cli::cli_div(theme = list(span.emph = list(color = "#ff4500")))
    cli::cli_abort("Comparison operator must be {.emph \"/\"} (e.g. {experiment}{.emph /}{control})")
  }
  check_samples(data, experiment, control)

  .term <- rlang::arg_match(.term, get_annotation_terms(data))
  .method <- rlang::arg_match(.method)
  if(is.null(data$analysis)){
    cli::cli_abort(c("No data for the analysis of expression differences.",
                     "i" = "  run expression() first."))
  }

  set_expression <- glue::glue("{experiment}/{control}")
  if(is.null(data$analysis[[set_expression]])){
    cli::cli_abort(c("No data for `{experiment}` `{control}`",
                     "i" = "  run expression({experiment}, {control}) first."))
  }

  data_expression <- data$analysis[[set_expression]]$expression

  data_expression_have <- names(data_expression)
  data_expression_need <- c(data$identifier, 'log2_foldchange')
  data_expression_diff <- setdiff(data_expression_need, data_expression_have)
  if(length(data_expression_diff) > 0) {
    cli::cli_div(theme = list(span.emph = list(color = "#ff4500"), span.info = list(color = "blue")))
    cli::cli_abort(c("x" = "Enrichment analysis data must contain the columns {.info {data_expression_need}}",
                     "i" = "Data is missing {.emph {data_expression_diff}}"))
  }

  if(.method == 'gsea') {
    table <- data_expression %>% enrichment_gsea(data, .term)
  } else if(.method == 'wilcoxon') {
    table <- data_expression %>% enrichment_wilcoxon(data, .term, ...)
  }

  data$analysis[[set_expression]]$enrichment[[.term]] <- table

  data$operations <- append(data$operations, glue::glue("Analysis: protein group enrichment via {.method}, grouping by {.term} for {experiment}/{control}"))

  return(data)
}
