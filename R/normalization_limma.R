#' Normalization function for a tidyproteomics data-object
#'
#' @param data tidyproteomics data object
#'
#' @return a tibble
#'
normalize_limma <- function(
    data = NULL
    ){

  # visible bindings
  identifier <- NULL

  nms <- c('identifier', 'sample', 'replicate', 'abundance_normalized')

    # create a data matrix for just, the values needed
  data_norm <- data %>%
    tidyr::pivot_wider(
      names_from = c("sample","replicate"),
      values_from = 'abundance',
      names_sep = "@"
    ) %>%
    as.data.frame() %>%
    tibble::column_to_rownames("identifier") %>%
    as.matrix() %>%
    limma::normalizeBetweenArrays(method="quantile") %>%
    as.data.frame() %>%
    tibble::rownames_to_column("identifier") %>%
    dplyr::as_tibble() %>%
    tidyr::pivot_longer(
      !identifier,
      names_to = 'sample',
      values_to = 'abundance_normalized'
    ) %>%
    tidyr::separate(sample, into = c("sample", "replicate"), sep="@") %>%
    dplyr::select(nms)

  return(data_norm)
}
