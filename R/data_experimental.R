#' Returns the data experimental set up
#'
#' @description
#' `experimental()` returns the transformative operations performed on the data.
#'
#' @param data tidyproteomics data object
#' @param destination a character string
#'
#' @return a character
#' @export
#'
#' @examples
#' library(dplyr, warn.conflicts = FALSE)
#' library(tidyproteomics)
#' #\dontrun{
#' hela_proteins <- path_to_package_data("p97KD_HCT116") %>%
#'    import("ProteomeDiscoverer", "proteins") %>%
#'    reassign(sample == "ctl", .replace = "control") %>%
#'    reassign(sample == "p97", .replace = "knockdown") %>%
#'    impute() %>%
#'    normalize(.method = c("linear","loess"))
#' }
#' hela_proteins %>% experimental()
#'
experimental <- function(
    data=NULL,
    destination=c('print','save')
){

  # visible bindings
  n <- NULL

  destination <- rlang::arg_match(destination)
  check_data(data)
  tbl <- data$experiments

  if(destination == "print") {
    cli::cli_h2(cli::style_bold("{.emph Experimental Design}"))
    print(as.data.frame(tbl), row.names = FALSE)
    println()
    glue::glue(" measurements:{nrow(tbl)} samples:{length(unique(tbl$sample))}")
  } else {
    tbl %>%
      readr::write_csv("table_data_experimental.csv")
  }
}

#' Main function for adding sample groups
#'
#' @param data a tidyproteomics data list-object
#' @param sample_groups a character string vector equal to the experimental row length
#'
#' @return a tidyproteomics data list-object
#' @export
#'

#'
experimental_groups <- function(
    data=NULL,
    sample_groups=NULL
){

  # visible bindings
  check_data(data)

  if(length(sample_groups) != nrow(data$experiments)){
    cli::cli_abort("sample_groups is length {length(sample_groups)}, and must be {nrow(data$experiments)}")
  }

  data$experiments <- data$experiments %>%
    dplyr::mutate(sample_group = sample_groups)

  return(data)
}
