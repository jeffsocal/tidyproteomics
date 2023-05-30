#' Main function for adding annotations to a tidyproteomics data-object
#'
#' @param data a tidyproteomics data list-object
#' @param annotations a character string vector
#' @param duplicates a character string, how to handle duplicate terms
#'
#' @return a tidyproteomics data list-object
#' @export
#'

#'
annotate <- function(
    data = NULL,
    annotations = NULL,
    duplicates = c("replace", "merge",  "leave")
){

  # visible bindings
  term <- NULL
  annotation <- NULL
  duplicate <- NULL

  check_data(data)
  duplicates <- rlang::arg_match(duplicates)
  if(is.null(annotations)) { return(data) }

  annotation_have <- colnames(annotations)
  annotation_need <- c(data$identifier, 'term', 'annotation')
  v_diff <- base::setdiff(annotation_need, annotation_have)

  if(data$analyte == 'peptides' & length(v_diff) > 0) {
    g_diff <- setdiff(data$identifier, v_diff)
    annotations <- annotations %>%
      dplyr::left_join(data$quantitative %>%
                  dplyr::select(data$identifier) %>%
                  unique(),
                by = g_diff)

    annotation_have <- colnames(annotations)
    annotation_need <- c(data$identifier, 'term', 'annotation')
    v_diff <- base::setdiff(annotation_need, annotation_have)
  }

  l_have <- annotation_have %>% length()
  l_need <- annotation_need %>% length()
  l_intr <- base::intersect(annotation_have, annotation_need) %>% length()

  if(l_have != l_need) {
    cli::cli_abort(c("x" = "Not all columns are present",
                     "i" = "Have {annotation_have}",
                     "i" = "Missing {v_diff}"))
  }

  data$annotations <- dplyr::bind_rows(
    data$annotations,
    annotations %>% dplyr::select(dplyr::all_of(annotation_need))) %>%
    unique() %>%
    dplyr::filter(!is.na(term)) %>%
    dplyr::filter(!is.na(annotation))

  if(duplicates == "merge") {
    data$annotations <- data$annotations %>%
      dplyr::group_by_at(c(data$identifier, 'term')) %>%
      dplyr::summarise(
        annotation = paste(annotation, collapse = "; "),
        .groups = 'drop'
      )
  } else if(duplicates == "replace") {
    data$annotations <- data$annotations %>%
      dplyr::group_by_at(c(data$identifier, 'term')) %>%
      dplyr::mutate(duplicate = dplyr::row_number()) %>%
      dplyr::slice_max(duplicate, n = 1, with_ties = FALSE) %>%
      dplyr::select(!duplicate) %>%
      dplyr::ungroup()
  }

  return(data)
}

#' Display the current annotation data
#'
#' @param data tidyproteomics data object
#' @param term a character string
#'
#' @return a vector
#' @export
#'
#' @examples
#' library(dplyr, warn.conflicts = FALSE)
#' library(tidyproteomics)
#'
#' hela_proteins %>% show_annotations()
#'
#' hela_proteins %>% show_annotations('reactome_pathway')
#'
show_annotations <- function(
    data,
    term = NULL
) {

  # visible bindings
  annotation <- NULL

  check_data(data)
  v_terms <- get_annotation_terms(data)

  if(is.null(term)) {
    cli::cli_inform(c("i" = "Current Annotation Terms"))
    cli::cli_ol()
    ulid <- cli::cli_ul()
    cli::cli_li(v_terms)
    cli::cli_end(ulid)
    cli::cli_end()
  } else {

    term <- rlang::arg_match(term, v_terms)
    anno <- get_annotations(data, term) %>%
      dplyr::select(annotation) %>%
      unique() %>%
      unlist()

    cli::cli_inform(c("i" = "Descriptions for {term}"))
    cli::cli_ol()
    ulid <- cli::cli_ul()
    cli::cli_li(anno)
    cli::cli_end(ulid)
    cli::cli_end()
  }

}
