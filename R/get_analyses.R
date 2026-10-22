#' Helper function to get all analyses
#'
#' @description
#' `get_analyses()` is a helper function that returns the summary or list of all
#' statistical analyses (differential expression, enrichment, etc.) stored within
#' a tidyproteomics data object.
#'
#' @param data tidyproteomics data object
#' @param format a character string specifying the return format: "table" (default) or "vector"
#'
#' @return a tibble or character vector
#' @export
#'
#' @examples
#' library(tidyproteomics)
#' hela_proteins %>%
#'   expression(knockdown/control) %>%
#'   get_analyses()
#'
get_analyses <- function(
    data = NULL,
    format = c("table", "vector")
) {
  check_data(data)
  format <- rlang::arg_match(format)

  if (is.null(data$analysis) || length(data$analysis) == 0) {
    if (format == "table") {
      return(tibble::tibble(
        comparison = character(),
        analysis = character(),
        term = character(),
        method = character()
      ))
    } else {
      return(character())
    }
  }

  rows <- list()
  for (comp in names(data$analysis)) {
    comp_obj <- data$analysis[[comp]]

    if (!is.null(comp_obj$expression)) {
      rows[[length(rows) + 1]] <- tibble::tibble(
        comparison = comp,
        analysis = "expression",
        term = NA_character_,
        method = "expression"
      )
    }

    if (!is.null(comp_obj$enrichment)) {
      for (term_name in names(comp_obj$enrichment)) {
        term_method <- comp_obj$enrichment[[term_name]]$method
        if (is.null(term_method)) {
          term_method <- "enrichment"
        }
        rows[[length(rows) + 1]] <- tibble::tibble(
          comparison = comp,
          analysis = "enrichment",
          term = term_name,
          method = as.character(term_method)
        )
      }
    }
  }

  if (length(rows) == 0) {
    tbl_out <- tibble::tibble(
      comparison = character(),
      analysis = character(),
      term = character(),
      method = character()
    )
  } else {
    tbl_out <- dplyr::bind_rows(rows)
  }

  if (format == "table") {
    return(tbl_out)
  } else {
    if (nrow(tbl_out) == 0) {
      return(character())
    }
    vec_out <- ifelse(
      is.na(tbl_out$term),
      paste(tbl_out$comparison, tbl_out$analysis, sep = ":"),
      paste(tbl_out$comparison, tbl_out$analysis, tbl_out$term, sep = ":")
    )
    return(as.character(vec_out))
  }
}
