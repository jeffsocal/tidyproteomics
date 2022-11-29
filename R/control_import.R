#' Main function for importing data to the tidyproteomics list object
#'
#' @param files a character vector of file paths
#' @param platform the source of the data (ProteomeDiscoverer, MaxQuant)
#' @param analyte the omics analyte (proteins, peptides)
#'
#' @return a tidyproteomics list data-object
#' @export
#'
#' @importFrom magrittr %>%
#'
#' @examples
#' library(dplyr, warn.conflicts = FALSE)
#' library(tidyproteomics)
#' ecoli_proteins <- path_to_package_data("ecoli-hint_proteins") %>%
#'    import("ProteomeDiscoverer", "proteins")
#'
#' ecoli_proteins %>% summary('sample')
#'
import <- function(
    files = NULL,
    platform = NULL,
    analyte = NULL
){

  # visible bindings
  data <- NULL

  if(is.null(files)) {cli::cli_abort(c("x" = "No files specified for importing."))}
  if(is.null(platform)) {cli::cli_abort(c("x" = "No file platform indicated"))}
  if(is.null(analyte)) {cli::cli_abort(c("x" = "No file analyte indicated"))}

  if(platform == 'ProteomeDiscoverer'){
    if(analyte == 'peptides') {data <- files %>% import_pd_peps()}
    if(analyte == 'proteins') {data <- files %>% import_pd_prots()}
    if(!analyte %in% c('peptides', 'proteins')) {cli::cli_abort(c("x" = "{analyte} not supported for {platform}"))}
  } else {
    data <- files %>% data_import(platform, analyte)
  }

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
