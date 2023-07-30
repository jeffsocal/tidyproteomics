#' helper function for having nice colors
#'
#' @return character vector of curated html colors
#'
theme_palette <- function(n = 16){

  pal <- c("#44709d", "#d97828", "#83992a", "#995d81",
    "#deb340", "#5c2751", "#a63d40",
    "#8ea4d2", "#704e2e", "#496f5d", "#e56399",
    "#f7b267", "#e4572e", "#69d1c5",
    "#48beff", "#8b8c89"  )

  if(n > 16){
    pal <- grDevices::rainbow(n = n, s = .67, v = .67)
  }

  return(pal)
}
