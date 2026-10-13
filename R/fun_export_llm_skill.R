#' Export the LLM agent skill document to the working directory
#'
#' @description
#' `export_llm_skill()` copies the bundled `SKILL.md` document from the package's
#' extdata directory to the current working directory or a specified destination.
#' This document contains complete instructions, API references, supported formats,
#' and workflows for use by Large Language Models (LLMs) and autonomous AI coding agents.
#'
#' @param destination path to directory where the file should be saved (default: ".")
#' @param file_name character string for the output filename (default: "SKILL.md")
#' @param overwrite boolean whether to overwrite an existing file (default: TRUE)
#'
#' @return invisible character string of the saved file path
#' @export
#'
#' @examples
#' \dontrun{
#' library(tidyproteomics)
#' export_llm_skill()
#' }
export_llm_skill <- function(
    destination = ".",
    file_name = "SKILL.md",
    overwrite = TRUE
){

  # Locate SKILL.md in installed package extdata
  source_path <- system.file("extdata", "SKILL.md", package = "tidyproteomics")

  # Fallback for development / uninstalled environments
  if (source_path == "" || !file.exists(source_path)) {
    candidates <- c(
      file.path("inst", "extdata", "SKILL.md"),
      file.path("..", "inst", "extdata", "SKILL.md")
    )
    for (candidate in candidates) {
      if (file.exists(candidate)) {
        source_path <- candidate
        break
      }
    }
  }

  if (source_path == "" || !file.exists(source_path)) {
    cli::cli_abort(c("x" = "SKILL.md file not found in package extdata"))
  }

  if (!dir.exists(destination)) {
    dir.create(destination, recursive = TRUE)
  }

  target_path <- file.path(destination, file_name)

  if (file.exists(target_path) && !overwrite) {
    cli::cli_abort(c("x" = "File {.file {target_path}} already exists and overwrite is FALSE."))
  }

  success <- file.copy(source_path, target_path, overwrite = overwrite)

  if (success) {
    cli::cli_alert_success("`{basename(target_path)}` exported to {destination}")
    invisible(normalizePath(target_path, mustWork = FALSE))
  } else {
    cli::cli_abort(c("x" = "Failed to copy SKILL.md to {.file {target_path}}"))
  }
}
