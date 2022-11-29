#' Helper function for saving plots
#'
#' @description
#' `plot_save` helper function
#'
#' @param plot a ggplot2 object
#' @param data a tidyproteomics data object
#' @param file_name a character string
#' @param destination a character string
#' @param height a numeric
#' @param width a numeric
#'
#' @return a ggplot2 object
#'
plot_save <- function(
    plot,
    data,
    file_name,
    destination = c("plot", "save", "png", "svg", "tiff", "jpeg"),
    height = 5,
    width = 8,
    ...
){

  # visible bindings
  ok_classes <- c('ggplot', 'eulergram', 'pheatmap')

  destination = rlang::arg_match(destination)
  if(!is.numeric(height) && (height < 1 | height > 30)) {cli::cli_abort("plot height must be a numeric between 1 and 30")}
  if(!is.numeric(width) && (width < 1 | width > 30)) {cli::cli_abort("plot width must be a numeric between 1 and 30")}
  if(length(intersect(class(plot), ok_classes)) == 0) {
    cli::cli_abort("plot must be one of type {ok_classes}")
  }

  if(destination == 'plot') {return(plot)}
  if(destination == 'save') {destination <- 'png'}

  file_name = glue::glue('plot_{file_name}.{destination}')

  plot %>% ggplot2::ggsave(filename=file_name, height = height, width = width, ...)
  cli::cli_alert_info("Saved {file_name}")
  return(data)
}
