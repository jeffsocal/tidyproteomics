#' Write table data locally
#'
#' @description
#' `save_table()` will save a summary tibble in the root directory of the
#' local project, based on the extension given in the file name. This is a
#' _passthrough_ function as it returns the original tibble.
#'
#' @param table a tibble
#' @param file_name a file name with extensions one of (.csv, .tsv, .rds, .xlsx)
#'
#' @return a tibble
#' @export
#'
#' @examples
#' library(dplyr, warn.conflicts = FALSE)
#' library(tidyproteomics)
#' hela_proteins %>%
#'    expression(knockdown/control) %>%
#'    export_analysis(knockdown/control, .analysis = "expression") %>%
#'    save_table("expression_limma_ko_over_wt.csv")
#'
save_table <- function(
    table,
    file_name = NULL
) {

  if(is.null(file_name)) {
    cli::cli_abort(c("x" = "No file name provided"))
  }

  format <- stringr::str_extract(file_name, "\\..{3,4}$")[1]
  base_name <- basename(file_name)
  path_name <- sub(base_name, "", file_name)

  format <- rlang::arg_match(format, c('.csv', '.xlsx','.rds'))
  cli::cli_process_start("Save table {base_name}")


  if(stringr::str_length(path_name) > 1 && !dir.exists(path_name)) {
    dir.create(path_name, recursive = T)
  }

  tryCatch({
    if(format == '.rds') {
      table %>% saveRDS(file_name)
    } else if(format == '.csv') {
      table %>% readr::write_csv(file_name)
    } else if(format == '.tsv') {
      table %>% readr::write_tsv(file_name)
    } else if(format == '.xlsx') {
      table %>% writexl::write_xlsx(file_name)
    }
  }, error = function(err) {
    err = as.character(as.vector(err))
    cli::cli_process_failed()
    cli::cli_abort(err)
  })

  cli::cli_process_done()
  return(table)
}

