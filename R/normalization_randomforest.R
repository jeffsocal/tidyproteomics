#' Normalization function for a tidyproteomics data-object
#'
#' @param data tidyproteomics list data-object
#' @param data_centered a tibble of centered values used for normalization
#'
#' @return a tibble
#'
normalize_randomforest <- function(
    data = NULL,
    data_centered = NULL
    ){

  # visible bindings
  identifier <- NULL
  abundance <- NULL
  abundance_centered <- NULL

  data_norm <- data %>%
    dplyr::group_by(sample, replicate) %>%
    tidyr::nest(data = c('identifier', 'abundance'))

  data_mdl <- data %>%
    dplyr::inner_join(data_centered, by='identifier', suffix = c("","_centered")) %>%
    dplyr::select(identifier, replicate, sample, abundance, abundance_centered) %>%
    dplyr::group_by(sample, replicate) %>%
    tidyr::nest(data = c('identifier', 'abundance', 'abundance_centered')) %>%
    dplyr::ungroup()

  for( i in 1:nrow(data_mdl) ){

    vals_medn <- data_mdl$data[[i]]$abundance_centered %>% unlist()
    vals_this <- data_mdl$data[[i]]$abundance %>% unlist()

    lm <- randomForest::randomForest(vals_medn ~ vals_this)

    data_norm$data[[i]]$abundance_normalized <- stats::predict(lm, newdata=data_norm$data[[i]]$abundance)
  }

  nms <- c('identifier', 'sample', 'replicate', 'abundance_normalized')

  data_norm <- data_norm %>%
    tidyr::unnest(data) %>%
    dplyr::ungroup() %>%
    dplyr::select(nms)

  return(data_norm)
}
