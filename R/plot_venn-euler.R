#' GGplot2 extension to plot a Euler diagram
#'
#' @param data a tidyproteomics data object
#' @param ... passthrough for ggsave see `plotting`
#'
#' @return a (tidyproteomics data-object | ggplot-object)
#' @export
#'
#' @examples
#' library(dplyr, warn.conflicts = FALSE)
#' library(tidyproteomics)
#' hela_proteins %>%
#'    subset(imputed == 0) %>%
#'    plot_euler()
#'
#' hela_proteins %>%
#'    subset(imputed == 0) %>%
#'    subset(cellular_component %like% "cytosol") %>%
#'    plot_euler()
#'
plot_euler <- function(
    data,
    ...){

  check_data(data)

  venn_set <- data %>% list_venn()

  plot <- plot(eulerr::euler(venn_set),
               main = glue::glue("Euler Diagram: {data$analyte}"),
               fills = list(fill = theme_palette(length(venn_set)), alpha = .33),
               quantities = T,
               legend = T,
               edges = T)

  return(plot_save(plot,
                   data,
                   glue::glue("{data$analyte}_euler"),
                   ...))
}

#' GGplot2 extension to plot a Venn diagram
#'
#' @param data a tidyproteomics data object
#' @param ... passthrough for ggsave see `plotting`
#'
#' @return a (tidyproteomics data-object | ggplot-object)
#' @export
#'
#' @examples
#' library(dplyr, warn.conflicts = FALSE)
#' library(tidyproteomics)
#' hela_proteins %>%
#'    subset(imputed == 0) %>%
#'    plot_venn()
#'
plot_venn <- function(
    data,
    ...){

  check_data(data)
  venn_set <- data %>% list_venn()

  plot <- plot(eulerr::venn(venn_set),
               main = glue::glue("Venn Diagram: {data$analyte}"),
               fills = list(fill = theme_palette(length(venn_set)), alpha = .33),
               quantities = T,
               legend = T,
               edges = T)

  return(plot_save(plot,
                   data,
                   glue::glue("{data$analyte}_venn"),
                   ...))

}

