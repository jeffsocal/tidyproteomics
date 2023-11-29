#' Helper function for iterative expression analysis
#'
#' @param data tidyproteomics data object
#' @param pairs the list of vector doublets
#'
#' @return list of vectors
#'
check_pairs <- function(
    pairs = NULL,
    sample_names = NULL
){

  if(is.null(sample_names)){ cli::cli_abort("No sample names given") }
  if(is.null(pairs)){ cli::cli_abort("No pairs given") }
  if(!is.list(pairs)){ cli::cli_abort("Pairs not supplied in 'list' format") }

  for(i in 1:length(pairs)){
    if(length(pairs[[i]]) != 2) { cli::cli_abort("Non-conforming pair at index {i}: {pairs[[i]]}") }
    if(length(intersect(sample_names, pairs[[i]])) != 2){
      for(n in 1:2){
        if(!pairs[[i]][n] %in% sample_names) cli::cli_abort("Unknown sample at index {i}: `{pairs[[i]][n]}`")
      }
    }
  }
}
