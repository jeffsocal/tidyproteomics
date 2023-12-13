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

#' Helper function extracting a subset of proteins
#'
#' @param data the tidyproteomics data object
#' @param include the set of proteins contained within the intersection of these samples
#' @param exclude the set of proteins found in these samples to exclude
#'
#' @return a character string
#'
intersect_venn <- function(
    data = NULL,
    include = NULL,
    exclude = NULL){

  check_data(data)
  sn <- data %>% get_sample_names()
  include <- rlang::arg_match(include, sn, multiple = T)
  if(!is.null(exclude))
    exclude <- rlang::arg_match(exclude, sn, multiple = T)

  venn_set <- data %>% list_venn()

  final <- venn_set[[include[1]]]
  if(length(include) > 1){
    for(i in 2:length(include)){
      final <- intersect(final, venn_set[[include[i]]])
    }
  }
  if(!is.null(exclude) & length(exclude) > 0){
    for(i in 1:length(exclude)){
      final <- setdiff(final, venn_set[[exclude[i]]])
    }
  }
  return(final)
}
