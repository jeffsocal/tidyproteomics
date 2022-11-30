#' reassign the sample info
#'
#' @description
#' `reassign()` enables editing of the sample descriptive in the experimental table.
#' This function will only replace the sample string and update the replicate number.
#'
#' @param data a tidyproteomics data-object
#' @param field a character string
#' @param pattern a character string
#' @param replace a character string
#'
#' @return a tidyproteomics data-object
#' @export
#'
#' @examples
#' library(dplyr, warn.conflicts = FALSE)
#' library(tidyproteomics)
#'
#' # check the experiment table
#' ecoli_proteins %>% summary("experiment")
#'
#' # make the modification
#' ecoli_proteins %>%
#'    reassign(field = "sample", pattern = "ko", replace = "hint_ko") %>%
#'    summary("sample")
#'
#' # reassign specific file_ids
#' ecoli_proteins %>%
#'    reassign(field = "sample_file", pattern = "f1", replace = "new") %>%
#'    reassign(field = "sample_file", pattern = "f2", replace = "new") %>%
#'    summary("sample")
#'
reassign <- function(
    data = NULL,
    field = NULL,
    pattern = NULL,
    replace = NULL
){

  # visible bindings
  sample_new <- NULL
  replicate_new <- NULL

  check_data(data)
  field <- rlang::arg_match(field, names(data$experiments))
  w <- which(grepl(pattern, unlist(data$experiments[,field])))

  if(grepl("^[0-9]", replace, perl = TRUE)){
    cli::cli_div(theme = list(span.emph = list(color = "#ff4500")))
    cli::cli_abort(c("x" = "Input {.emph sample} names must not start with a {.emph numeric}"))
  }

  if(grepl("/", replace, perl = TRUE)){
    cli::cli_div(theme = list(span.emph = list(color = "#ff4500")))
    cli::cli_abort(c("x" = "Input {.emph sample} names must not contain a {.emph /}"))
  }

  if(is.null(replace)){
    if(length(w) == 0) {
      cli::cli_abort("did not find {.emph pattern} in {.emph field}`")
    } else {
      cli::cli_h2("found the following {length(w)} rows")
      print(as.data.frame(data$experiments[w,]), row.names = FALSE)
    }
  } else if(!is.character(replace)){
    cli::cli_abort("`replace` must be a charcter string")
  } else if(length(w) == 0) {
    cli::cli_abort("did not find {.emph {pattern}} in {.emph {field}}")
  }else {

    w <- which(grepl(pattern, as.data.frame(data$experiments)[,field]))
    if(length(w) == 0) {return(data)}

    data$experiments$sample[w] <- replace
    data$experiments <- data$experiments %>%
      dplyr::group_by(sample) %>%
      dplyr::mutate(replicate = as.character(dplyr::row_number())) %>%
      dplyr::ungroup()

    data$quantitative <- data$quantitative %>%
      dplyr::select(!c("sample","replicate")) %>%
      dplyr::full_join(
        data$experiments %>% dplyr::select("sample_id","sample","replicate"),
        by = "sample_id"
      ) %>%
      dplyr::relocate(c('sample_id', 'sample', 'replicate', data$identifier))

    data$operations <- append(data$operations, glue::glue("Data reassigned sample to {replace} where '{field} like {pattern}'"))
    return(data)
  }
}
