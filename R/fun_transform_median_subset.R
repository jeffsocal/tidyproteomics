#' helper function for normalizing quantitative data from a tidyproteomics data-object
#'
#' @param data tidyproteomics data object
#' @param group_by character vector
#' @param rename character string
#'
#' @importFrom rlang :=
#'
#' @return a tibble
#'
transform_median <- function(
    data,
    group_by = c('identifier'),
    rename = "log2_med"
){

  # visible bindings
  abundance <- NULL

  data  %>%
    dplyr::group_by_at(group_by) %>%
    dplyr::summarise(
      abundance = stats::median(abundance, na.rm = T),
      .groups = "drop"
    ) %>%
    dplyr::filter(!is.na(abundance)) %>%
    dplyr::rename({{rename}} := abundance)
}


#' helper function for normalizing quantitative data from a tidyproteomics data-object
#'
#' @param data tidyproteomics data object
#' @param data_factor tidyproteomics data object
#' @param ... pass through arguments
#'
#' @return a tibble
#'
transform_factor <- function(
    data,
    data_factor = NULL,
    ...
){

  data_med <- data %>% transform_log2() %>% transform_median(...)

  if(!is.null(data_factor)){
    pre_n <- data_med$identifier %>% unique() %>% length()
    pre_range <- data_med$log2_med %>% range()

    data_med <- data_factor %>% transform_log2() %>% transform_median(...)

    pst_n <- data_med$identifier %>% unique() %>% length()
    pst_range <- data_med$log2_med %>% range()

    cli::cli_alert_warning(" normalization based on {pst_n} of {pre_n} identifiers")
    if(pre_range[1] < pst_range[1] | pre_range[2] > pst_range[2]){
      cli::cli_alert_warning("  {.emph WARNING}: filter narrowed range, NAs may result")
    }
  }
  return(data_med)
}
