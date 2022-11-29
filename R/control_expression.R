#' Summarize the data
#'
#' @description
#' `expression()` is an analysis function that computes the protein summary
#' statistics for a given tidyproteomics data object.
#'
#' @param data tidyproteomics data object
#' @param ... two sample comparison e.g. experimental/control
#' @param .method a two-distribution test function returning a p_value for the null
#' hypothesis. Example functions include t.test, wilcox.test, stats::ks.test,
#' additionally, the string _"limma"_ can be used to select from the limma
#' package to compute an empirical Bayesian estimation which performs better with
#' non-linear distributions and uneven replicate balance between samples.
#'
#' @return a tibble
#' @export
#'
#' @examples
#' library(dplyr, warn.conflicts = FALSE)
#' library(tidyproteomics)
#'
#' # simple t.test expression analysis
#' ecoli_proteins %>%
#'    expression(ko/wt) %>%
#'    export_analysis(ko/wt, .analysis = "expression")
#'
#' # a wilcox.test expression analysis
#' ecoli_proteins %>%
#'    expression(ko/wt, .method = stats::wilcox.test) %>%
#'    export_analysis(ko/wt, .analysis = "expression")
#'
#' # a one-tailed wilcox.test expression analysis
#' wilcoxon_less <- function(x, y) {
#'    stats::wilcox.test(x, y, alternative = 'less')
#' }
#' ecoli_proteins <- ecoli_proteins %>%
#'    expression(ko/wt, .method = stats::wilcox.test)
#'
#' ecoli_proteins %>% export_analysis(ko/wt, .analysis = "expression")
#'
#' # Note: the userdefined function is preserved in the operations tracking
#' ecoli_proteins %>% operations()
#'
#' # limma expression analysis
#' ecoli_proteins %>%
#'    expression(ko/wt, .method = "limma") %>%
#'    export_analysis(ko/wt, .analysis = "expression")
#'
expression <- function(
    data = NULL,
    ...,
    .method = stats::t.test
){

  # visible bindings
  imputed <- NULL
  log2_foldchange <- NULL
  foldchange <- NULL

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

  if(mode(.method) == 'function') {
    table <- data %>% expression_test(experiment, control, .method)
    .method <- gsub("\\.", "_", as.character(methods::functionBody(.method))[2])
  } else if(mode(.method) == 'character' && .method == 'limma') {
    table <- data %>% expression_limma(experiment, control)
  } else {
    cli::cli_div(theme = list(span.emph = list(color = "#ff4500")))
    cli::cli_abort(c("x" = "invalid test function {.emph {.method}}"))
  }

  sample_ids = data$experiment %>%
    dplyr::filter(sample %in% c(experiment, control)) %>%
    dplyr::select(sample_id) %>%
    unlist()

  # add in the counts and imputation stats
  if(length(base::intersect(c('match_between_runs', 'imputed'), names(data$accounting))) > 0){
    table <- data$accounting %>%
      dplyr::filter(sample_id %in% sample_ids) %>%
      tidyr::pivot_longer(dplyr::matches("match|impute"), names_to = 'type', values_to = 'imputed') %>%
      dplyr::filter(!is.na(imputed)) %>%
      dplyr::group_by(dplyr::across(c(data$identifier, "sample_id"))) %>%
      dplyr::summarise(imputed = min(imputed) == 1, .groups = 'drop') %>%
      dplyr::group_by(dplyr::across(data$identifier)) %>%
      dplyr::summarise(imputed = sum(imputed)/dplyr::n(), n = dplyr::n()) %>%
      dplyr::left_join(table, by = data$identifier)
  }

  set_expression <- glue::glue("{experiment}/{control}")
  data$analysis[[set_expression]]$expression <- table %>%
    dplyr::arrange(dplyr::desc(log2_foldchange)) %>%
    dplyr::filter(!is.na(foldchange))

  data$operations <- append(data$operations, glue::glue("Analysis: expression difference {.method} {experiment}/{control}"))

  return(data)
}

check_samples <- function(
    data,
    experiment = NULL,
    control = NULL
){
  experiment <- rlang::arg_match(experiment, unique(data$experiments$sample))
  control <- rlang::arg_match(control, unique(data$experiments$sample))

  if(control == experiment){
    cli::cli_div(theme = list(span.emph = list(color = "#ff4500"), span.info = list(color = "blue")))
    cli::cli_abort("x" = "Expression analysis must have different samples choosen",
                   "{.info experiment}:{.emph {experiment}} and {.info control}:{.emph {control}} are the same")
  }
}
