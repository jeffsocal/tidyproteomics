#' Helper function to get all accounting terms
#'
#' @param data tidyproteomics data object
#'
#' @return a vector
#'
get_accountings <- function(
    data = NULL
) {

  names_accounting <- names(data$accounting)
  names_accounting <- names_accounting[-which(grepl("sample|protein$|peptide$|modification$|group$", names_accounting))]

  return(names_accounting)
}
