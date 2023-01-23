#' Main method for imputing missing values
#'
#' @param data a tidyproteomics list data-object
#' @param impute_function summary statistic function. Default is base::min, examples of
#' other functions include min, max, mean, sum. One could write a function for the
#' lower 5%, lower5th <- function(x) { quantile(x, 0.05)[[1]] }. Note, NAs will be
#' be removed in the function call.
#' @param minimum_to_impute the minimum ratio to impute at
#' @param cores the number of threads used to speed the calculation
#'
#' @return a tidyproteomics list data-object
#' @export
#'
#' @examples
#' library(dplyr, warn.conflicts = FALSE)
#' library(tidyproteomics)
#' hela_proteins %>% summary("sample")
#'
#' hela_proteins %>% impute(impute_function = stats::median) %>% summary("sample")
#'
#' hela_proteins %>% impute(impute_function = "randomforest") %>% summary("sample")
#'
impute <- function(
    data = NULL,
    impute_function = base::min,
    method = c('within', 'between'),
    minimum_to_impute = 0.25,
    cores = 2
){

  # visible bindings
  imputed <- NULL
  imputed.x <- NULL
  imputed.y <- NULL

  check_data(data)
  method <- rlang::arg_match(method)
  quant_source <- data$quantitative_source
  table <- data %>% extract(quant_source, na.rm = TRUE)

  if('imputed' %in% colnames(data$accounting)) {
    cli::cli_alert_info("This is data has been collapsed from the peptide-level")
    cli::cli_alert_info("... imputation should be done at that level")
  }

  cli::cli_progress_bar(type = 'tasks')

  if(mode(impute_function) == 'function') {

    impute_function_str <- rlang::quo_text(impute_function)
    cli::cli_div(theme = list(span.emph = list(color = "#ff4500")))
    cli::cli_progress_step("Imputing {.emph {method}} samples by {.emph {impute_function_str}}")

    table <- table %>% impute_byfunction(impute_function, minimum_to_impute = minimum_to_impute, method = method)
    data$operations <- append(data$operations, glue::glue("Missing values imputed {method} samples via {impute_function_str}."))

  } else if(mode(impute_function) == 'character' && impute_function == 'randomforest') {

    method <- 'between'

    cli::cli_div(theme = list(span.emph = list(color = "#ff4500")))
    cli::cli_progress_step("Imputing {.emph {method}} samples by {.emph randomforest}")

    l_out <- table %>% impute_randomforest(minimum_to_impute = minimum_to_impute, cores = cores)
    table <- l_out$table
    data$operations <- append(data$operations, glue::glue("Missing values imputed {method} samples via {impute_function}."))
    data$operations <- append(data$operations, l_out$operation)
  } else {
    cli::cli_div(theme = list(span.emph = list(color = "#ff4500")))
    cli::cli_abort(c("x" = "invalid impute_function {.emph {impute_function}}"))
  }

  cli::cli_progress_done()

  table <- data$experiments %>%
    dplyr::select(c('sample_id', 'sample', 'replicate')) %>%
    dplyr::full_join(table, by = c('sample', 'replicate'))

  data <- data %>%
    merge_quantitative(table %>% dplyr::select(!imputed), quant_source)

  table <- table %>%
    munge_identifier("separate", data$identifier) %>%
    dplyr::select(!dplyr::matches('abundance')) %>%
    dplyr::select(!dplyr::matches('sample_id')) %>%
    dplyr::full_join(data$experiments, by = c('sample', 'replicate')) %>%
    dplyr::select(!c('sample', 'replicate', 'import_file', 'sample_file'))

  join_names <- colnames(table)[-which("imputed" == colnames(table))]
  data$accounting <- data$accounting %>%
    dplyr::full_join(table, by = join_names)

  if('imputed.y' %in% colnames(data$accounting)){
    data$accounting <- data$accounting %>%
      dplyr::rename(imputed = imputed.y) %>%
      dplyr::select(!c('imputed.x'))
  }

  impute_msg <- glue::glue("... {data$accounting %>% dplyr::filter(imputed == TRUE) %>% nrow()} values imputed")
  cli::cli_alert_info(impute_msg)
  data$operations <- append(data$operations, impute_msg)

  return(data)
}
