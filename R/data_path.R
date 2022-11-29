#' Helper function for displaying path to data
#'
#' @param item a character string
#'
#' @return print the table to console
#' @export
#'
path_to_package_data <- function(
    item = c('proteins','peptides', 'fasta')
){

  item <- item[1]

  configs <- list(
    list.files(system.file("extdata/config", "", package = "tidyproteomics"), full.names = T, pattern = item),
    list.files(system.file("extdata", "", package = "tidyproteomics"), full.names = T, pattern = item)
  ) %>% unlist()

  if(is.na(configs[1])) {
    cli::cli_div(theme = list(span.emph = list(color = "#ff4500")))
    cli::cli_abort(c("x" = "No path found for {.emph {item}}"))
  }

  sub("\\/{2,}", "/", configs[1])
}
