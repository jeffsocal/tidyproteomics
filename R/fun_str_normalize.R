#' Normalize the column names in a tibble
#'
#' @param x a vector
#'
#' @return a vector
#' @export
#'
str_normalize <- function(x) {
  x <- base::tolower(x)
  x <- gsub("\\:|\\[|\\]|\\(|\\)|\\,|\\.|\\/", " ", x)
  x <- gsub("\\#", "num", x)
  x <- gsub("\\%", "percent", x)
  x <- trimws(x)
  x <- gsub("\\s+", "_", x)
  return(x)
}
