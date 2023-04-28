#' Normalization function for a tidyproteomics data-object
#'
#' @param data tidyproteomics list data-object
#' @param data_centered a tibble of centered values used for normalization
#'
#' @return a tibble
#'
normalize_scaled <- function(
    data = NULL,
    data_centered = NULL
){

  # visible bindings
  abundance_all <- NULL
  abundance <- NULL

  data_centered <- data_centered %>%
    dplyr::left_join(data %>%
                        center(group_by = c('sample', 'replicate'), values = 'abundance', method='sum') %>%
                        dplyr::rename(abundance = dplyr::matches('abundance')),
                      by = c('sample', 'replicate'), suffix = c('','_all')) %>%
    dplyr::mutate(ratio = abundance_all / abundance)

  t_median <- data_centered$abundance %>% stats::median(na.rm = T)

  data_norm <- list()
  for( i in 1:nrow(data_centered) ){
    tdf <- data %>% dplyr::filter(sample == data_centered$sample[i] & replicate == data_centered$replicate[i])
    t_factor <- t_median * data_centered$ratio[i]
    data_norm[[i]] <- tdf %>%
      dplyr::mutate(abundance_normalized = (abundance / sum(abundance, na.rm = T)) * t_factor)
  }

  nms <- c('identifier', 'sample', 'replicate', 'abundance_normalized')

  data_norm <- data_norm %>%
    dplyr::bind_rows() %>%
    dplyr::select(nms)

  return(data_norm)
}
