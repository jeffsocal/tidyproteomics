#' Helper function to export the config file to current project directory
#'
#' @param platform the source of the data (ProteomeDiscoverer, MaxQuant)
#' @param analyte the omics analyte (proteins, peptides)
#'
#' @return success or fail
#' @export
#'
#' @importFrom magrittr %>%
#'
#' @examples
#' library(tidyproteomics)
#' #\dontrun{
#' export_config("mzTab", 'peptides')
#' }
#'
export_config <- function(
    platform = NULL,
    analyte = c('proteins','peptides')
){

  # visible bindings
  if(is.null(platform)) {cli::cli_abort(c("x" = "No file platform indicated"))}

  analyte <- rlang::arg_match(analyte)

  l_config <- list.files(system.file("extdata/config", "",
                                     package = "tidyproteomics"), full.names = F)
  platform <- rlang::arg_match(platform, sub("\\_(proteins|peptides)\\.tsv", "", l_config) %>% unique())

  path_to_config <- glue::glue("{platform}_{analyte}") %>% path_to_package_data()

  if(!file.exists(path_to_config)) {cli::cli_abort("Config file not found for {platform}_{analyte}")}

  file.copy(path_to_config, basename(path_to_config))
  cli::cli_alert_success("`{basename(path_to_config)}` exported to project directory")
}
