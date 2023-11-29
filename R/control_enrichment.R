#' Compute protein enrichment
#'
#' @description
#' `enrichment()` is an analysis function that computes the protein summary
#' statistics for a given tidyproteomics data object.
#'
#' @param data tidyproteomics data object
#' @param ... two sample comparison e.g. experimental/control
#' @param .pairs a list of vectors each containing two named sample groups
#' @param .term a character string referencing ".term" in the annotations table
#' @param .method a character string
#' @param .score_type a character string. From the fgsea manual: "This parameter
#' defines the GSEA score type. Possible options are ("std", "pos", "neg"). By
#' default ("std") the enrichment score is computed as in the original GSEA. The
#' "pos" and "neg" score types are intended to be used for one-tailed tests
#' (i.e. when one is interested only in positive ("pos") or negateive ("neg") enrichment)."
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
#'    enrichment(knockdown/control, .term = "biological_process") %>%
#'    export_analysis(knockdown/control, .analysis = "enrichment", .term = "biological_process")
#'
#' # using a Wilcoxon Rank Sum method
#' hela_proteins %>%
#'    expression(knockdown/control) %>%
#'    enrichment(knockdown/control, .term = "biological_process", .method = "wilcoxon") %>%
#'    export_analysis(knockdown/control, .analysis = "enrichment", .term = "biological_process")
#'
#' # using the .pairs argument when multiple comparisons are needed
#' comps <- list(c("control","knockdown"),
#'             c("knockdown","control"))
#'
#' hela_proteins %>%
#'    expression(.pairs = comps) %>%
#'    enrichment(.pairs = comps, .term = "biological_process")
#'
enrichment <- function(
    data = NULL,
    ...,
    .pairs = NULL,
    .term = NULL,
    .method = c('gsea', 'wilcoxon'),
    .score_type = c("std", "pos", "neg"),
    .cpu_cores = 1
){

  check_data(data)
  str_quo <- tidyproteomics_quo(...)
  if(!is.null(.pairs)) {
    cli::cli_inform("Using the supplied {length(.pairs)} sample pairs ...")
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
  cli::cli_progress_bar(type = 'tasks')

  for(i in 1:length(.pairs)){

    ui_t <- .term <- rlang::arg_match(.term, get_annotation_terms(data))
    ui_m <- .method <- rlang::arg_match(.method)
    ui_s <- .score_type <- rlang::arg_match(.score_type)

    experiment <- .pairs[[i]][1]
    control <- .pairs[[i]][2]
    cli::cli_div(theme = list(span.emph = list(color = "#ff4500")))
    cli::cli_progress_step(" .. enrichment::{.emph {ui_m}} testing {.emph {experiment} / {control}} by term {.emph {ui_t}}")

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

    if(.method == 'gsea') {
      table <- data_expression %>% enrichment_gsea(data, .term, .score_type, .cpu_cores)
      data$operations <- append(data$operations, glue::glue("Analysis: protein group enrichment via {.method}::{.score_type}, grouping by {.term} for {experiment}/{control}"))
    } else if(.method == 'wilcoxon') {
      table <- data_expression %>% enrichment_wilcoxon(data, .term, cpu_cores = .cpu_cores)
      data$operations <- append(data$operations, glue::glue("Analysis: protein group enrichment via {.method}, grouping by {.term} for {experiment}/{control}"))
    }

    data$analysis[[set_expression]]$enrichment[[.term]] <- table

    cli::cli_progress_done()
  }

  return(data)
}
