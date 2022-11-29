#' Check the integrity of a tidyproteomics data object
#'
#' @description
#' `check_data()` is a helper function that checks the structure and contents of
#' a tidyproteomics data object
#'
#' @param data tidyproteomics data object
#'
#' @return silent on success, an abort message on fail
#'
check_data <- function(
    data=NULL
){

  # fail if data is NULL
  if(is.null(data)) {
    cli::cli_div(theme = list(span.emph = list(color = "#ff4500")))
    cli::cli_abort(c("x" = "Input cannot be {.emph NULL}"))
    }
  if(class(data) != 'tidyproteomics') {
    cli::cli_div(theme = list(span.emph = list(color = "#ff4500")))
    cli::cli_abort(c("x" = "Input must be of type {.emph tidyproteomics}"))
  }

  data_elements_need <- c("origin","analyte","identifier","quantitative_source","operations",
                          "experiments","annotations","quantitative","accounting","analysis")
  data_elements <- names(data)

  # fail if analyte is not a named element
  if(!"analyte" %in% data_elements) {
    cli::cli_div(theme = list(span.emph = list(color = "#ff4500")))
    cli::cli_abort(c("x" = "Input {.emph data-object} not recognized"))
  }
  analyte <- data$analyte
  # analyte <- rlang::arg_match(analyte, c('peptides','proteins','sequences'))

  if(data$analyte == 'sequences'){
    sequence_elements_need <- c("protein","sequence","coverage","residues","peptides","modifications")
    sequence_elements <- names(data$quantitative[[1]])
    sequence_elements <- rlang::arg_match(sequence_elements, sequence_elements_need, multiple = T)
  } else {
    data_elements <- rlang::arg_match(data_elements, data_elements_need, multiple = TRUE)

    if(length(which(grepl("^[0-9]", data$experiments$sample, perl = TRUE))) > 0){
      cli::cli_div(theme = list(span.emph = list(color = "#ff4500")))
      cli::cli_abort(c("x" = "Input {.emph sample} names must not start with a {.emph numeric}"))
    }

    if(length(which(grepl("\\-", data$experiments$sample, perl = TRUE))) > 0){
      cli::cli_div(theme = list(span.emph = list(color = "#ff4500")))
      cli::cli_abort(c("x" = "Input {.emph sample} names must not contain a {.emph hyphen (-)}"))
    }
  }
}

#' Create a crc32 hash on a vector
#'
#' @description
#' `hash_vector()` is a helper function that returns a crc32 hash on a vector
#'
#' @param x a vector
#'
#' @return a hash of x
#'
hash_vector <- function(x){ unlist(lapply(x, digest::digest, 'crc32'))}
