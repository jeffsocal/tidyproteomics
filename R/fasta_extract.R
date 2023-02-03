#' Get the string defined by the regex
#'
#' @description
#' `fasta_extract()` get the current string based on regex
#'
#' @param string a character
#' @param regex a list
#'
#' @return a list
#'
fasta_extract <- function(
    string = NULL,
    regex = NULL
){
  string <- unlist(strsplit(string, "\n"))
  out <- list()
  for( i in 1:length(regex) ) {
    out_name <- names(regex[i])
    if(out_name == 'sequence'){
      out_str <- paste(string[-1], collapse = "")
    } else {
      out_str <- unlist(stringr::str_extract(string[1], regex[[i]]))
    }
    out[[out_name]] <- out_str
  }
  return(out)
}
