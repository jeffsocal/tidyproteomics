#' Helper function to get all annotations for a given term
#'
#' @param data tidyproteomics data object
#' @param term a character string
#' @param format a character string specifying the return format: "table" (default) or "vector"
#'
#' @return a tibble or character vector
#' @export
#'
#' @examples
#' library(tidyproteomics)
#' get_annotations(hela_proteins, "biological_process")
#' get_annotations(hela_proteins, "biological_process", format = "vector")
#'
get_annotations <- function(
    data = NULL,
    term = NULL,
    format = c("table", "vector")
) {

  # visible bindings
  term_group <- term
  annotation <- NULL

  check_data(data)
  format <- rlang::arg_match(format)

  tbl <- data$annotations %>%
    dplyr::filter(term == term_group) %>%
    tidyr::separate_rows(annotation, sep = ";") %>%
    dplyr::mutate(annotation = trimws(annotation)) %>%
    dplyr::select(!term)

  if (format == "table") {
    return(tbl)
  }

  return(unique(tbl$annotation))
}

#' Helper function to get available terms
#'
#' @param data tidyproteomics data object
#' @param format a character string specifying the return format: "table" (default) or "vector"
#'
#' @return a tibble or character vector
#' @export
#'
#' @examples
#' library(tidyproteomics)
#' get_annotation_terms(hela_proteins)
#' get_annotation_terms(hela_proteins, format = "vector")
#'
get_annotation_terms <- function(
    data = NULL,
    format = c("table", "vector")
) {
  check_data(data)
  format <- rlang::arg_match(format)

  vec_terms <- c(data$identifier, unique(data$annotations$term))

  if (format == "table") {
    id_name <- data$identifier[1]
    id_count <- length(unique(unlist(data$quantitative[[id_name]])))

    rows <- list(
      tibble::tibble(
        term = id_name,
        type = "identifier",
        unique_annotations = id_count
      )
    )

    if (!is.null(data$annotations) && nrow(data$annotations) > 0) {
      anno_terms_ordered <- unique(data$annotations$term)
      anno_counts <- data$annotations %>%
        dplyr::group_by(term) %>%
        dplyr::summarise(
          unique_annotations = dplyr::n_distinct(annotation),
          .groups = "drop"
        ) %>%
        dplyr::arrange(match(term, anno_terms_ordered)) %>%
        dplyr::mutate(type = "annotation") %>%
        dplyr::relocate(term, type, unique_annotations)

      rows[[2]] <- anno_counts
    }

    return(dplyr::bind_rows(rows))
  }

  return(vec_terms)
}

