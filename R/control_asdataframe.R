#' Helper function to convert the data-object into a tibble
#'
#' @description
#' `as.data.frame()` is a function that converts the tidyproteomics data object into
#' a tibble. This tibble is in the long-format, such that a there is a single
#' observation per line.
#'
#' @param data tidyproteomics data object
#' @param shape the orientation of the quantitative data as either a single measure
#' per row (**long**), or as multiple measures per protein/peptide (**wide**).
#' @param values indicates the selected normalization to output. The default is
#' that selected at the time of normalization.
#' @return a tibble
#' @exportS3Method
#'
#' @examples
#' library(dplyr, warn.conflicts = FALSE)
#' library(tidyproteomics)
#'
#' # convert the data-object to a data.frame
#' hela_proteins %>% as.data.frame() %>% as_tibble()
#'
#' # select the wide format
#' hela_proteins %>% as.data.frame(shape = 'wide') %>% as_tibble()
#'
#' # select the wide format & drop some columns
#' hela_proteins %>%
#'    as.data.frame(shape = 'wide',
#'                  drop = c('description','wiki_pathway','reactome_pathway','biological_process')) %>%
#'    as_tibble()
#'
as.data.frame.tidyproteomics <- function(
    data,
    shape = c('long', 'wide'),
    values = NULL,
    drop = NULL
){

  check_data(data)
  if(is.null(values)) { values <- data$quantitative_source }
  if(!is.null(drop)) {
    can_drop <- c(colnames(data$accounting), data$annotations$term |> unique())
    if(length(intersect(drop, can_drop)) == 0){
      cli::cli_abort("trying to drop `{drop}`, only finding `{can_drop}`")
    }
  }
  values <- rlang::arg_match(values, get_quant_names(data))
  shape <- rlang::arg_match(shape)

  # pull the data from each table
  tbl_exp <- data$experiments %>% dplyr::select(!"import_file")
  tbl_acc <- data$accounting
  tbl_ann <- data$annotations
  tbl_qnt <- data$quantitative %>%
    dplyr::select(!tidyselect::matches('abundance') | tidyselect::matches(paste0(values, "$")))

  if(!is.null(tbl_ann)) {
    tbl_ann <- tbl_ann %>% tidyr::pivot_wider(names_from = 'term', values_from = 'annotation')
  }

  # build the output table
  tbl_out <- tbl_exp %>%
    dplyr::inner_join(tbl_qnt, by = c('sample_id', 'sample', 'replicate'))

  if(!is.null(tbl_ann)) { tbl_out <- tbl_out %>% dplyr::left_join(tbl_ann, by = data$identifier) }
  if(!is.null(tbl_acc)) { tbl_out <- tbl_out %>% dplyr::left_join(tbl_acc, by = c(data$identifier, 'sample_id')) }

  # pivot to wide if requested
  if(shape == 'wide'){

    pivot_on_values <- paste0(values, "$")

    if(!'match_between_runs' %in% drop & 'match_between_runs' %in% colnames(data$accounting)) {
      pivot_on_values <- paste0(pivot_on_values, "|match_between_runs")
    }
    if(!'imputed' %in% drop & 'imputed' %in% colnames(data$accounting)) {
      pivot_on_values <- paste0(pivot_on_values, "|imputed")
    }

    drop <- c('sample_id', 'sample_file', drop)

    tbl_out <- tbl_out %>%
      dplyr::select(!dplyr::all_of(drop)) %>%
      dplyr::mutate(origin = values) %>%
      tidyr::pivot_wider(
        names_from = c('sample', 'replicate', 'origin'),
        values_from = tidyselect::matches(pivot_on_values)
      )
  } else if(!is.null(drop)) {
    tbl_out <- tbl_out %>% dplyr::select(!dplyr::all_of(drop))
  }

  return(base::as.data.frame(tbl_out))
}
