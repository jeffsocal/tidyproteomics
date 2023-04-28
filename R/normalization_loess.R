#' Normalization function for a tidyproteomics data-object
#'
#' @param data tidyproteomics list data-object
#' @param data_centered a tibble of centered values used for normalization
#'
#' @return a tibble
#'
normalize_loess <- function(
    data = NULL,
    data_centered = NULL
){

  # visible bindings
  identifier <- NULL
  replicate <- NULL
  sample <- NULL
  abundance <- NULL
  abundance_centered <- NULL

  pre_range <- data$abundance %>% range()
  pst_range <- data_centered$abundance %>% range(na.rm = TRUE)

  adjust_loess_control <- F
  if(pre_range[1] < pst_range[1] | pre_range[2] > pst_range[2]){
    adjust_loess_control <- T
  }

  data_norm <- data %>%
    dplyr::left_join(data_centered, by='identifier', suffix = c("","_centered")) %>%
    dplyr::select(identifier, replicate, sample, abundance, abundance_centered) %>%
    dplyr::group_by(sample, replicate) %>%
    tidyr::nest(data = c('identifier', 'abundance', 'abundance_centered')) %>%
    dplyr::mutate(m = purrr::map(data, stats::loess, formula = abundance_centered ~ abundance))

  for( i in 1:nrow(data_norm) ){
    tdf <- data_norm$data[[i]]
    lm <- data_norm$m[[i]]
    if(adjust_loess_control == T){
      data_norm$data[[i]]$abundance_normalized <- stats::predict(lm, newdata=tdf$abundance, control=stats::loess.control(surface="direct"))
    } else {
      data_norm$data[[i]]$abundance_normalized <- stats::predict(lm, newdata=tdf$abundance)
    }
  }

  nms <- c('identifier', 'sample', 'replicate', 'abundance_normalized')

  data_norm <- data_norm %>%
    tidyr::unnest(data) %>%
    dplyr::ungroup() %>%
    dplyr::select(nms)

  return(data_norm)
}
