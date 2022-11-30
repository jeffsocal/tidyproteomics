#' Plot normalized values
#'
#' @description
#' `plot_normalization()` is a GGplot2 implementation for plotting the normalization
#' effects visualized as a box plot.
#'
#' @param data tidyproteomics data object
#' @param ... passthrough for ggsave see `plotting`
#'
#' @return a (tidyproteomics data-object | ggplot-object)
#' @export
#'
#' @examples
#' library(dplyr, warn.conflicts = FALSE)
#' library(tidyproteomics)
#' ecoli_proteins %>%
#'   normalize(.method = c("scaled", "median", "linear", "limma", "loess")) %>%
#'   plot_normalization()
#'
plot_normalization <- function(
    data = NULL,
    ...
){

  # visible bindings
  origin <- NULL
  abundance <- NULL

  check_data(data)

  analyte <- data$analyte
  quant_values_all <- c('raw', 'median', 'scaled', 'linear',
                        'limma', 'loess', 'svm', 'randomforest')
  quant_values <- get_quant_names(data)
  data_quant <- data %>% extract(values = quant_values)

  norm_vals <- intersect(quant_values_all, quant_values)

  data_medians <- data_quant %>%
    dplyr::group_by(origin) %>%
    dplyr::summarise(abundance = stats::  median(abundance, na.rm=T), .groups = 'drop')

  plot <- data_quant %>%
    dplyr::mutate(origin = forcats::fct_relevel(origin, norm_vals)) %>%
    ggplot2::ggplot() +
    ggplot2::geom_hline(data = data_medians, ggplot2::aes(yintercept = abundance)) +
    ggplot2::geom_boxplot(ggplot2::aes(sample, abundance, color=replicate)) +
    ggplot2::scale_y_log10() +
    # ggplot2::scale_color_manual(values = theme_palette) +
    ggplot2::facet_wrap(~origin, nrow=1) +
    ggplot2::coord_flip() +
    ggplot2::theme_minimal() +
    ggplot2::scale_color_manual(values = theme_palette()) +
    ggplot2::labs(title = "Normalization", subtitle = data$analyte)

  if(length(norm_vals) > 3){
    plot <- plot +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 30, hjust=1, size=8))
  }

  return(plot_save(plot,
                   data,
                   glue::glue("{data$analyte}_normalization"),
                   ...))
}
