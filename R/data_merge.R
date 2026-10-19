#' Merge two or more tidyproteomics data objects
#'
#' @description
#' `merge()` returns a single tidyproteomics data object from multiple.
#'
#' @param x a tidyproteomics data object
#' @param y a tidyproteomics data object or list of tidyproteomics data objects
#' @param ... additional tidyproteomics data objects to merge
#' @param quantitative_source a character string indicating which quantitative
#' value to merge on ('raw', 'selected', or 'all'). If `selected` is chosen then each
#' dataset's specific normalization will be used and renamed to 'abundance_selected'.
#' If `all` is chosen, then the possibility exists that some normalization values
#' will fill in with NAs.
#'
#' @return a tidyproteomics data object
#' @exportS3Method base::merge
#'
merge.tidyproteomics <- function(
    x,
    y = NULL,
    ...,
    quantitative_source = c('raw', 'selected', 'all')
){
  quantitative_source <- match.arg(quantitative_source)
  data_list <- list(x)
  if (!is.null(y)) {
    if (inherits(y, "tidyproteomics")) {
      data_list <- c(data_list, list(y))
    } else if (is.list(y)) {
      data_list <- c(data_list, y)
    }
  }
  dots <- list(...)
  for (item in dots) {
    if (inherits(item, "tidyproteomics")) {
      data_list <- c(data_list, list(item))
    } else if (is.list(item)) {
      data_list <- c(data_list, item)
    }
  }
  merge_tidyproteomics(data_list = data_list, quantitative_source = quantitative_source)
}

#' Merge multiple tidyproteomics data-objects
#'
#' @description
#' `merge_tidyproteomics()` returns a single tidyproteomics data object from multiple.
#'
#' @param data_list a list of tidyproteomics data objects
#' @param quantitative_source a character string indicating which quantitative
#' value to merge on. If `selected` is chosen then each dataset's specific normalization
#' will be used and renamed to 'abundance_selected'. If `all` is chosen, then the
#' possibility exists that some normalization values will fillin with NAs.
#'
#' @return a tidyproteomics data object
#' @export
#'
merge_tidyproteomics <- function(
    data_list = NULL,
    quantitative_source = c('raw', 'selected', 'all')
){

  if (inherits(data_list, "tidyproteomics")) {
    data_list <- list(data_list)
  }

  if (!is.list(data_list) || length(data_list) < 2) {
    cli::cli_abort("merge_tidyproteomics requires a list of at least 2 tidyproteomics objects")
  }

  quantitative_source <- rlang::arg_match(quantitative_source)

  l_operations <- list()
  l_experiments <- list()
  l_quantitative <- list()
  l_accouting <- list()
  l_annotations <- list()

  for(i in 1:length(data_list)){

    this_data <- data_list[[i]]
    check_data(this_data)
    if(i == 1) {
      this_analyte <- this_data$analyte
      this_identifier <- this_data$identifier
    } else {
      if(this_analyte != this_data$analyte) {
        cli::cli_abort(c("x" = "Analytes of the first set `{this_analyte}` do not match the second `{this_data$analyte}`"))
      }
    }
    this_quantitative_source <- quantitative_source
    if(quantitative_source != 'all'){

      quant_source <- 'raw'
      if(quantitative_source == 'selected') {
        this_quantitative_source <- 'selected'
        quant_source <- this_data$quantitative_source
      }

      this_data$quantitative <- this_data$quantitative %>%
        dplyr::select(!tidyselect::matches('abundance') | tidyselect::matches(paste0(quant_source, "$")))

      if(quantitative_source == 'selected') {
        this_data$quantitative <- this_data$quantitative %>%
          dplyr::rename(abundance_selected = tidyselect::matches(quant_source))
      }
    }

    l_operations <- append(l_operations, paste0(this_data$origin, " [", i, "]: ", this_data$operations))
    l_experiments[[i]] <- this_data$experiments
    l_quantitative[[i]] <- this_data$quantitative
    l_accouting[[i]] <- this_data$accounting
    l_annotations[[i]] <- this_data$annotations

  }

  new_data <- list()
  new_data$origin <- "Merged"
  new_data$analyte <- this_analyte
  new_data$identifier <- this_identifier
  new_data$quantitative_source <- this_quantitative_source
  new_data$operations <- append(l_operations, glue::glue("Merged {length(data_list)} data sets"))
  new_data$experiments <- dplyr::bind_rows(l_experiments)
  new_data$quantitative <- dplyr::bind_rows(l_quantitative)
  new_data$accounting <- dplyr::bind_rows(l_accouting)

  tbl_annotations <- dplyr::bind_rows(l_annotations) %>% unique()
  if(nrow(tbl_annotations) == 0) tbl_annotations <- NULL
  new_data$annotations <- tbl_annotations

  class(new_data) <- 'tidyproteomics'
  return(new_data)
}


#' Helper function merging normalized data back into the main data-object
#'
#' @param data tidyproteomics data subset tibble
#' @param data_quant tidyproteomics data subset tibble
#' @param values character string vector
#'
#' @return a tibble
#'
merge_quantitative <- function(
    data = NULL,
    data_quant = NULL,
    values = 'raw'
){

  check_data(data)

  data_quant <- data_quant %>% munge_identifier("separate", data$identifier)

  colnames(data_quant)[which('abundance' == names(data_quant))] <- paste0('abundance_', values)

  # fixed issue with matrix based computations reintroducing NA's into the
  # quantitative data table
  data$quantitative <- data$quantitative %>%
    dplyr::select(!dplyr::matches(paste0("\\_", values))) %>%
    dplyr::full_join(
      data_quant, by = names(data_quant)[which(!grepl('abundance', names(data_quant)))]
    ) %>%
    dplyr::filter(!is.na(sample_id))

  return(data)
}
