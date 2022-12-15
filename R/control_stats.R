#' Summarize the data
#'
#' @description
#' `summary()` is an analysis function that computes the protein summary
#' statistics for a given tidyproteomics data object. This is a _passthrough_ function
#' as it returns the original tidyproteomics data-object.
#'
#' @param data tidyproteomics data object
#' @param by what to summarize
#' @param destination character string, one of (save, print)
#' @param contamination as character string
#'
#' @return a tibble on *print*, a tidyproteomics data-object on *save*
#' @export
#'
#' @examples
#' library(dplyr, warn.conflicts = FALSE)
#' library(tidyproteomics)
#'
#' # a global summary
#' hela_proteins %>% summary()
#'
#' # a summary by sample
#' hela_proteins %>% summary("sample")
#'
#' # a summary by sample with match_between_runs removed
#' hela_proteins %>%
#'    subset(match_between_runs == FALSE) %>%
#'    summary("sample")
#'
#' # a summary of match_between_runs
#' hela_proteins %>% summary("match_between_runs")
#'
#' hela_proteins %>% summary("cellular_component")
#'
#' hela_proteins %>% summary("biological_process")
#'
summary.tidyproteomics <- function(
    data,
    by = c('global'),
    destination = c("print", "save", "return"),
    limit = 25,
    contamination = NULL
){

  check_data(data)
  destination <- rlang::arg_match(destination)
  if(!is.numeric(limit)) {cli::cli_abort("limit must be a numeric not `{limit}`")}
  limit <- max(ceiling(limit), 1)

  if(is.null(by)) {by <- 'global'}

  if(!is.null(contamination)) {
    by <- 'contamination'
    table <- data %>% stats_contamination(pattern = contamination)
  } else if(by == 'global') {
    table <- data %>% stats_summary()
  } else {

    vars <- data %>% get_unique_variables(by)
    if(is.null(vars)) {return(NULL)}

    if(is.numeric(vars) == TRUE) {
      vars <- sort(vars, decreasing = TRUE)
    } else {
      vars <- sort(vars)
    }

    if(length(vars) > limit) {
      cli::cli_alert_info("Too many variables, limiting to the first {limit}")
      vars <- vars[1:limit]
    }

    table <- list()
    for( var in vars ) {
      table[[rlang::quo_text(var)]] <- data %>%
        subset(!!by == !!var, .verbose = FALSE, rm.mbr = FALSE) %>%
        stats_summary() %>%
        dplyr::mutate(variable = var)
    }
    table <- table %>% dplyr::bind_rows() %>%
      dplyr::relocate(variable)
    colnames(table)[1] <- by
  }

  if(destination == 'save'){
    table %>% readr::write_csv(paste0("table_", data$analyte, "_summary_", by,".csv"))
    return(data)
  } else if(destination == 'return'){
    return(table)
  } else {
    stats_print(table, paste("Summary:", by))
  }
}
