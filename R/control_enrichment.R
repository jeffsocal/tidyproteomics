#' Compute protein enrichment
#'
#' @description
#' `enrichment()` is an analysis function that computes the protein summary
#' statistics for a given tidyproteomics data object.
#'
#' @param data tidyproteomics data object
#' @param ... two sample comparison e.g. experimental/control
#' @param .pairs a list of vectors each containing two named sample groups
#' @param .terms a character string referencing "term(s)" in the annotations table
#' @param .method a character string
#' @param .score_type a character string. From the fgsea manual: "This parameter
#' defines the GSEA score type. Possible options are ("std", "pos", "neg"). By
#' default ("std") the enrichment score is computed as in the original GSEA. The
#' "pos" and "neg" score types are intended to be used for one-tailed tests
#' (i.e. when one is interested only in positive ("pos") or negateive ("neg") enrichment)."
#' @param .log2fc_min used only for Fisher's Exact Test, a numeric defining the minimum log2 foldchange to consider as "enriched"
#' @param .significance_max used only for Fisher's Exact Test, a numeric defining the maximum statistical significance to consider as "enriched"
#' @param .cpu_cores the number of threads used to speed the calculation
#'
#' @return a tibble
#' @export
#'
#' @examples
#' library(dplyr, warn.conflicts = FALSE)
#' library(tidyproteomics)
#'
#' # using the default GSEA method
#' hela_proteins %>%
#'    expression(knockdown/control) %>%
#'    enrichment(knockdown/control, .terms = "biological_process") %>%
#'    export_analysis(knockdown/control, .analysis = "enrichment", .term = "biological_process")
#'
#' # using a Wilcoxon Rank Sum method
#' hela_proteins %>%
#'    expression(knockdown/control) %>%
#'    enrichment(knockdown/control, .terms = "biological_process", .method = "wilcoxon") %>%
#'    export_analysis(knockdown/control, .analysis = "enrichment", .term = "biological_process")
#'
#' # using the .pairs argument when multiple comparisons are needed
#' comps <- list(c("control","knockdown"),
#'             c("knockdown","control"))
#'
#' hela_proteins %>%
#'    expression(.pairs = comps) %>%
#'    enrichment(.pairs = comps, .terms = c("biological_process", "molecular_function")
#'
enrichment <- function(
    data = NULL,
    ...,
    .pairs = NULL,
    .terms = NULL,
    .method = c('gsea', 'wilcoxon', 'fishers_exact'),
    .score_type = c("std", "pos", "neg"),
    .log2fc_min = 0,
    .significance_min = 0.05,
    .cpu_cores = 1
){

  check_data(data)
  str_quo <- tidyproteomics_quo(...)
  if(!is.null(.pairs)) {
    cli::cli_inform("Enrichment Analysis - using the supplied {length(.pairs)} sample pairs ...")
  } else if(is.null(str_quo)) {
    return(data)
  } else {
    .pairs <- list(c(str_quo[['variable']], str_quo[['value']]))
    if(str_quo[['operator']] != "/") {
      cli::cli_div(theme = list(span.emph = list(color = "#ff4500")))
      cli::cli_abort("Comparison operator must be {.emph \"/\"} (e.g. {experiment}{.emph /}{control})")
    }
  }

  # a quick check on pairs to mitigate any issues prior to computing
  check_pairs(.pairs, get_sample_names(data))
  identifier <- data$identifier
  avaiable_terms <- get_annotation_terms(data)
  cli::cli_progress_bar(type = 'tasks')

  for(i in 1:length(.pairs)){

    ui_m <- .method <- rlang::arg_match(.method)
    ui_s <- .score_type <- rlang::arg_match(.score_type)

    experiment <- .pairs[[i]][1]
    control <- .pairs[[i]][2]
    cli::cli_div(theme = list(span.emph = list(color = "#ff4500")))

    check_samples(data, experiment, control)

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


    for( use_term in .terms){

      if(!use_term %in% avaiable_terms) { next() }

      cli::cli_progress_step(" .. enrichment::{.emph {ui_m}} testing {.emph {experiment} / {control}} by term {.emph {use_term}}")

      if(.method == 'gsea') {
        table <- data_expression %>% enrichment_gsea(data, use_term, .score_type, .cpu_cores)
        str_method = glue::glue("GSEA::{.score_type}")
      } else if(.method == 'wilcoxon') {
        table <- data_expression %>% enrichment_wilcoxon(data, use_term, cpu_cores = .cpu_cores)
        str_method = "Wilcoxon"
      } else if(.method == 'fishers_exact') {
        table <- data_expression %>% enrichment_fishersexact(data, use_term,
                                                             log2fc_min = .log2fc_min,
                                                             significance_min = .significance_min,
                                                             cpu_cores = .cpu_cores)
        str_method = "Fisher's Exact"
      }
        data$operations <- append(data$operations, glue::glue("Analysis: protein group enrichment via {str_method}, grouping by {use_term} for {experiment}/{control}"))

      data$analysis[[set_expression]]$enrichment[[use_term]] <- list(method = glue::glue("{str_method} grouping by {use_term}"),
                                                                     data = table)
      cli::cli_progress_done()

    }

  }

  return(data)
}
