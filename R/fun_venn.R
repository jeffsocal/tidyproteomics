#' Helper function Venn and Euler plots
#'
#' @param data tidyproteomics data object
#' @param ... pass through arguments
#'
#' @return list of vectors
#'
list_venn <- function(
    data = NULL,
    ...
){

  # visible bindings
  identifier <- NULL
  abundance <- NULL

  # use just to test the data
  check_data(data)
  data_quant <- data %>%
    extract() %>%
    dplyr::group_by(identifier, sample) %>%
    dplyr::slice_max(abundance, n=1, with_ties = FALSE) %>%
    dplyr::select(sample, identifier) %>%
    dplyr::group_by(sample) %>%
    tidyr::nest()

  vset <- list()
  # for(i in 1:nrow(data_quant)){
  #   vset[[i]] <- data_quant$data[[i]] %>%
  #     unlist() %>% as.character()
  # }

  for(samplei in data_quant$sample){
    vset[[samplei]] <- data_quant %>%
      dplyr::filter(sample == samplei) %>%
      dplyr::ungroup() %>%
      dplyr::select(data) %>%
      unlist() %>% as.character()
  }

  return(vset)
}
