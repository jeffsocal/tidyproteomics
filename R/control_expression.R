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
#' @param .p.adjust a stats::p.adjust string for multiple test correction, default is
#' 'BH' (Benjamini & Hochberg, 1995)
#'
#' @return a tibble
#' @export
#'
#' @examples
#' library(dplyr, warn.conflicts = FALSE)
#' library(tidyproteomics)
#'
#' # simple t.test expression analysis
#' hela_proteins %>%
#'    expression(knockdown/control) %>%
#'    export_analysis(knockdown/control, .analysis = "expression")
#'
#' # a wilcox.test expression analysis
#' hela_proteins %>%
#'    expression(knockdown/control, .method = stats::wilcox.test) %>%
#'    export_analysis(knockdown/control, .analysis = "expression")
#'
#' # a one-tailed wilcox.test expression analysis
#' wilcoxon_less <- function(x, y) {
#'    stats::wilcox.test(x, y, alternative = "less")
#' }
#' hela_proteins <- hela_proteins %>%
#'    expression(knockdown/control, .method = stats::wilcox.test)
#'
#' hela_proteins %>% export_analysis(knockdown/control, .analysis = "expression")
#'
#' # Note: the userdefined function is preserved in the operations tracking
#' hela_proteins %>% operations()
#'
#' # limma expression analysis
#' hela_proteins %>%
#'    expression(knockdown/control, .method = "limma") %>%
#'    export_analysis(knockdown/control, .analysis = "expression")
#'
#' # using the .pairs argument when multiple comparisons are needed
#' comps <- list(c("control","knockdown"),
#'             c("knockdown","control"))
#'
#' hela_proteins %>%
#'    expression(.pairs = comps)
#'
expression <- function(
    data = NULL,
    ...,
    .pairs = NULL,
    .method = stats::t.test,
    .p.adjust = 'BH'
){

  # visible bindings
  imputed <- NULL
  log2_foldchange <- NULL
  foldchange <- NULL
  sample_id <- NULL

  check_data(data)
  str_quo <- tidyproteomics_quo(...)
  if(!is.null(.pairs)) {
    cli::cli_inform("Expression Analysis - using the supplied {length(.pairs)} sample pairs ...")
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
  cli::cli_progress_bar(type = 'tasks')

  for(i in 1:length(.pairs)){

    experiment <- .pairs[[i]][1]
    control <- .pairs[[i]][2]

    check_samples(data, experiment, control)

    if(mode(.method) == 'function') {
      this_method <- gsub("\\.", "_", as.character(methods::functionBody(.method))[2])
      cli::cli_div(theme = list(span.emph = list(color = "#ff4500")))
      cli::cli_progress_step(" .. expression::{.emph {this_method}} testing {.emph {experiment} / {control}}")

      table <- data %>% expression_test(experiment, control, .method = .method, .p.adjust = .p.adjust)
    } else if(mode(.method) == 'character' && .method == 'limma') {
      this_method <- .method
      cli::cli_div(theme = list(span.emph = list(color = "#ff4500")))
      cli::cli_progress_step(" .. expression::{.emph {this_method}} testing {.emph {experiment} / {control}}")

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

      table_imput_class <- data$accounting %>%
        dplyr::filter(sample_id %in% sample_ids) %>%
        dplyr::left_join(data$experiment |> dplyr::select(sample_id, sample), by = "sample_id") %>%
        tidyr::pivot_longer(dplyr::matches("match|impute"), names_to = 'type', values_to = 'imputed') %>%
        dplyr::filter(!is.na(imputed)) %>%
        dplyr::group_by_at(c(data$identifier, "sample_id", "sample")) %>%
        dplyr::summarise(imputed = min(imputed) == 1, .groups = 'drop') %>%
        dplyr::group_by_at(c(data$identifier, 'sample')) %>%
        dplyr::summarise(imputed = sum(imputed)/dplyr::n()) %>%
        tidyr::pivot_wider(
          names_from = 'sample', values_from = 'imputed', names_prefix = "imputed_"
        )

      table <- data$accounting %>%
        dplyr::filter(sample_id %in% sample_ids) %>%
        tidyr::pivot_longer(dplyr::matches("match|impute"), names_to = 'type', values_to = 'imputed') %>%
        dplyr::filter(!is.na(imputed)) %>%
        dplyr::group_by_at(c(data$identifier, "sample_id")) %>%
        dplyr::summarise(imputed = min(imputed) == 1, .groups = 'drop') %>%
        dplyr::group_by_at(data$identifier) %>%
        dplyr::summarise(imputed = sum(imputed)/dplyr::n(),
                         n = dplyr::n(), .groups = 'drop') %>%
        dplyr::left_join(table_imput_class, by = data$identifier) %>%
        dplyr::left_join(table, by = data$identifier)
    }

    set_expression <- glue::glue("{experiment}/{control}")
    data$analysis[[set_expression]]$expression <- table %>%
      dplyr::arrange(dplyr::desc(log2_foldchange)) %>%
      dplyr::filter(!is.na(foldchange))

    data$operations <- append(data$operations, glue::glue("Analysis: expression difference {this_method} {experiment}/{control}, p.adjust = {.p.adjust}"))

    cli::cli_progress_done()
  }
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
