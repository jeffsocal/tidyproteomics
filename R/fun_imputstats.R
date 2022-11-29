#' Helper function for calculating imputation stats
#'
#' @param x a tibble
#'
#' @return list of vectors
#'
impute_ratio <- function(x){
  #visible bindings
  imputed <- NULL
  w <- x %>%
    dplyr::group_by(sample) %>%
    dplyr::summarise(
      ratio = length(which(imputed == TRUE)) / length(imputed),
      .groups = 'drop'
    )
  return(max(w$ratio))
}
