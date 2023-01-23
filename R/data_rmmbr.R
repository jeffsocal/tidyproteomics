#' Remove MBR from the dataset across segments
#'
#' @description
#' `rm.mbr()` function is designed to remove match_between_runs between segments.
#' This function will return a smaller tidyproteomics data-object.
#'
#' @param data tidyproteomics data object
#' @param ... a three part expression (eg. x == a)
#' @param .groups a character string
#'
#' @return a tibble
#' @export
#'
#' @examples
#' library(dplyr, warn.conflicts = FALSE)
#' library(tidyproteomics)
#'
#' hela_proteins %>%
#'    summary('sample')
#'
#' hela_proteins %>%
#'    rm.mbr(.groups = 'sample') %>%
#'    summary('sample')
#'
rm.mbr <- function(
    data = NULL,
    ...,
    .groups = c('all', 'sample')
){

  # visible bindings
  match_between_runs <- NULL
  n_mbr <- NULL

  check_data(data)
  .groups <- rlang::arg_match(.groups)

  n_ids_pre <- data$quantitative %>%
    dplyr::select(data$identifier) %>%
    unique() %>% nrow()

  v_groups <- data$identifier
  if(.groups == "sample")
    v_groups <- c(v_groups, 'sample')

  r2a <- data$accounting %>%
    dplyr::inner_join(data$experiments, by = 'sample_id') %>%
    dplyr::group_by(dplyr::across(v_groups)) %>%
    dplyr::summarise(
      n_mbr = length(which(match_between_runs == FALSE)),
      .groups = 'drop'
    ) %>%
    dplyr::filter(n_mbr != 0) %>%
    dplyr::select(v_groups)

  data$quantitative <- data$quantitative %>%
    dplyr::inner_join(r2a, by = v_groups)

  data$accounting <- data$quantitative %>%
    dplyr::select(c("sample_id", data$identifier)) %>%
    dplyr::inner_join(data$accounting, by = c("sample_id", data$identifier))

  n_ids_pst <- data$quantitative %>%
    dplyr::select(data$identifier) %>%
    unique() %>% nrow()

  data$operations <- append(data$operations, glue::glue("Removed MBR between `{.groups}` groups"))

  return(data)
}
