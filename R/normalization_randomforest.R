#' Normalization function for a tidyproteomics data-object
#'
#' @param data tidyproteomics list data-object
#' @param data_centered a tibble of centered values used for normalization
#'
#' @return a tibble
#'
normalize_randomforest <- function(
    data = NULL,
    data_centered = NULL,
    .cores = 1
){

  # visible bindings
  identifier <- NULL
  abundance <- NULL
  abundance_centered <- NULL

  data_norm <- data %>%
    dplyr::group_by(sample, replicate) %>%
    tidyr::nest(data = c('identifier', 'abundance'))

  data_centered <- data_centered %>%
    dplyr::filter(!is.na(abundance))

  data_mdl <- data %>%
    dplyr::inner_join(data_centered, by='identifier', suffix = c("","_centered")) %>%
    dplyr::select(identifier, replicate, sample, abundance, abundance_centered) %>%
    dplyr::group_by(sample, replicate) %>%
    tidyr::nest(data = c('identifier', 'abundance', 'abundance_centered')) %>%
    dplyr::ungroup()


  # if only one core is selected don't bother with the parallel package
  if(.cores > 1){
    data_mdl$data <- parallel::mclapply(data_mdl$data,
                                        rf_parallel,
                                        mc.cores = .cores)
  } else {
    data_mdl$data <- lapply(data_mdl$data,
                            rf_parallel)
  }

  # old for loop
  # for( i in 1:nrow(data_mdl) ){
  #
  #   vals_medn <- data_mdl$data[[i]]$abundance_centered %>% unlist()
  #   vals_this <- data_mdl$data[[i]]$abundance %>% unlist()
  #
  #   nlm <- randomForest::randomForest(vals_medn ~ vals_this)
  #
  #   data_norm$data[[i]]$abundance_normalized <- stats::predict(nlm, newdata=data_norm$data[[i]]$abundance)
  # }

  nms <- c('identifier', 'sample', 'replicate', 'abundance_normalized')

  data_norm <- data_mdl %>%
    tidyr::unnest(data) %>%
    dplyr::ungroup() %>%
    dplyr::select(nms)

  return(data_norm)
}

#' parallel compute function for randomforest
#'
#' @param df a tibble of raw and centered values
#'
#' @return a tibble
#'
rf_parallel <- function(df){
  vals_medn <- df$abundance_centered %>% unlist()
  vals_this <- df$abundance %>% unlist()

  nlm <- randomForest::randomForest(vals_medn ~ vals_this)

  df$abundance_normalized <- stats::predict(nlm, newdata=vals_this)

  return(df)
}
