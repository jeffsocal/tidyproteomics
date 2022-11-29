#' Helper function for displaying data
#'
#' @param table a tibble
#' @param title a character string
#'
#' @return print the table to console
#'
stats_print <- function(
    table,
    title = NULL
){

  if(!tibble::is_tibble(table)) { return() }

  cli::cli_end()
  if(!is.null(title)) {
    cli::cli_h2(title)
  }
  table %>% as.data.frame() %>% print(row.names = FALSE)
  cli::cli_text("")
}
