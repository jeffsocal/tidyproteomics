#' set a named vector
#'
#' @param config a data.frame of configuration values
#' @param category a character string
#'
#' @return a named vector
#'
set_vect <- function(
    config = NULL,
    category = NULL
){
  if(is.null(category)) {
    return(stats::setNames(config$column_import, config$column_defined))
  }

  w <- which(config$category == category)
  if(length(w) == 0) {return(NULL)}
  stats::setNames(config$column_import[w], config$column_defined[w])
}

#' match a named vector to string vector
#'
#' @param un_vec an un-named vector
#' @param n_vec a named vector
#'
#' @return a named vector
#'
match_vect <- function(
    un_vec,
    n_vec
){

  g_vec <- c()
  for(n_int in n_vec){
    g_vec <-   n_vec
  }
  g_vec <- intersect(n_vec, un_vec)
  # names(g_vec) <- names(n_vec[which(n_vec %in% g_vec)])
  return(g_vec)
}
