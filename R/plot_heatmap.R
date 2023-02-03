#' Plot a heatmap of quantitative values by sample
#'
#' @description
#' `plot_heatmap()` is a pheatmap implementation for plotting the commonly
#' visualized quantitative heatmap according to sample. Both the samples and the
#' quantitative values are clustered and visualized.
#'
#' @param data tidyproteomics data object
#' @param tag a character string
#' @param row_names a boolean
#' @param ... passthrough for ggsave see `plotting`
#'
#' @return a (tidyproteomics data-object | ggplot-object)
#' @export
#'
#' @examples
#' library(dplyr, warn.conflicts = FALSE)
#' library(tidyproteomics)
#' hela_proteins %>%
#'   normalize(.method = c("scaled", "median", "linear", "limma", "loess")) %>%
#'   select_normalization() %>%
#'   plot_heatmap()
#'
plot_heatmap <- function(
    data = NULL,
    tag=NULL,
    row_names=FALSE,
    ...
){

  # visible bindings
  color <- NULL
  abundance <- NULL
  sample_rep <- NULL

  check_data(data)

  quantval <- data$quantitative_source
  data_quant <- data %>% extract(values = quantval, na.rm = TRUE) %>%
    dplyr::select(!dplyr::matches('^origin$'))
  theme_palette <- theme_palette()

  data_col <- data_quant %>%
    dplyr::group_by(sample) %>%
    dplyr::summarise(
      n = replicate %>% unique() %>% length(),
      .groups = 'drop'
    )

  data_col <- data_col %>%
    dplyr::mutate(color = theme_palette[1:nrow(data_col)])

  group_color <- list(sample = data_col %>% dplyr::select(color) %>% unlist() %>% as.vector())
  names(group_color$sample) <- data_col$sample

  data_munge <- data_quant %>%
    dplyr::mutate(abundance = log10(abundance)) %>%
    tidyr::pivot_wider(
      names_from = c('sample', 'replicate'),
      values_from = 'abundance',
      values_fill = 0
    ) %>% as.data.frame() %>%
    tibble::column_to_rownames('identifier') %>%
    as.matrix()

  col <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(10, "RdYlBu"))(32) %>% rev()

  plot <- pheatmap::pheatmap(data_munge,
                             scale = "none",
                             silent = TRUE,
                             col =  col,
                             show_rownames = row_names,
                             annotation_colors = group_color,
                             #annotation_legend = FALSE,
                             #legend = FALSE,
                             annotation_names_col = FALSE,
                             annotation_col = data_quant %>%
                               dplyr::select(sample, replicate) %>%
                               unique() %>%
                               tidyr::unite(sample_rep, sample, replicate, remove = F) %>%
                               dplyr::select(!replicate) %>%
                               as.data.frame() %>%
                               tibble::column_to_rownames('sample_rep'),
                             main = glue::glue("Heatmap: {data$analyte} {data$quantitative_source} \n")
  )

  if(!'null device' %in% names(grDevices::dev.cur())) { invisible(grDevices::dev.off()) }

  return(plot_save(plot,
                   data,
                   glue::glue("{data$analyte}_{data$quantitative_source}_heatmap"),
                   ...))
}
