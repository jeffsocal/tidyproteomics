#' Merge multiple tidyproteomics data-objects
#'
#' @description
#' `merge()` returns a single tidyproteomics data object from multiple.
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
merge <- function(
    data_list = NULL,
    quantitative_source = c('raw', 'selected', 'all')
){

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
    this_quantitative_source <- 'raw'
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
  new_data$annotations <- dplyr::bind_rows(l_annotations) %>% unique()

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

  data$quantitative <- data$quantitative %>%
    dplyr::select(!dplyr::matches(paste0("\\_", values))) %>%
    dplyr::full_join(
      data_quant, by = names(data_quant)[which(!grepl('abundance', names(data_quant)))]
    )

  return(data)
}
