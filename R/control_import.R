#' Main function for importing data
#'
#' @description `import()` reads files from various platforms into the
#' tidyproteomics data object -- see also the documentation `vignette("importing")`
#' and `vignette("workflow-importing")`
#'
#' @param files a character vector of file paths
#' @param platform the source of the data (ProteomeDiscoverer, MaxQuant, etc.)
#' @param analyte the omics analyte (proteins, peptides)
#' @param path a character string pointing to the local configuration file (directory/file.tsv)
#'
#' @return a tidyproteomics list data-object
#' @export
#'
#' @importFrom magrittr %>%
#'
#' @examples
#' library(dplyr, warn.conflicts = FALSE)
#' library(tidyproteomics)
#' hela_proteins <- path_to_package_data("p97KD_HCT116") %>%
#'    import("ProteomeDiscoverer", "proteins")
#' hela_proteins %>% summary("sample")
#'

import <- function(
    files = NULL,
    platform = NULL,
    analyte = NULL,
    path = NULL
){

  # visible bindings
  data <- NULL

  if(is.null(files)) {cli::cli_abort(c("x" = "No files specified for importing."))}
  if(is.null(platform)) {cli::cli_abort(c("x" = "No file platform indicated"))}
  if(is.null(analyte)) {cli::cli_abort(c("x" = "No file analyte indicated"))}

  data <- files %>% data_import(platform = platform, analyte = analyte, path = path)

  if(nrow(data$quantitative) == 0) {
    cli::cli_abort(c("x" = "... not data present"))
  }

  data$operations <- list(glue::glue("Data files ({basename(files)}) were imported as {analyte} from {platform}"))

  return(tidyproteomics(data))
}

#' Tidy-Quant data object print definition
#'
#' @param obj tidyproteomics data object
#'
#' @return print object summary
#'

#' @export
#'
tidyproteomics <- function(obj) {
  class(obj) <- "tidyproteomics"
  return(obj)
}
