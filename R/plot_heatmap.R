#' Plot a heatmap of quantitative values by sample
#'
#' @description
#' `plot_heatmap()` is a pheatmap implementation for plotting the commonly
#' visualized quantitative heatmap according to sample. Both the samples and the
#' quantitative values are clustered and visualized.
#'
#' @param data tidyproteomics data object
#' @param tag a character string
#' @param row_names either FALSE, or an annotation (eg. protein, gene_name, description etc.)
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
  row_names_bool <- FALSE

  check_data(data)

  quantval <- data$quantitative_source
  data_quant <- data %>%
    extract(values = quantval, na.rm = TRUE) %>%
    dplyr::select(!dplyr::matches('^origin$'))


  if(is.character(row_names)){
    row_names_bool <- TRUE
    pos_names <- get_annotation_terms(data)
    if(!row_names %in% pos_names){
      cli::cli_abort("`row_names` must be one of {pos_names}")
    }
  } else {
    row_names_bool <- row_names
    row_names <- data$identifier[1]
  }

  if( row_names != data$identifier ){
    tbl_annotation <- data |> get_annotations(row_names)
    id_col <- which(colnames(tbl_annotation) %in% data$identifier)
    colnames(tbl_annotation)[id_col] <- 'identifier'
    id_col <- which(colnames(tbl_annotation) == 'annotation')
    colnames(tbl_annotation)[id_col] <- row_names

    data_quant <- data_quant |>
      dplyr::inner_join(
        tbl_annotation, by = 'identifier'
      )

    id_col <- which(colnames(data_quant) == 'identifier')
    data_quant <- data_quant[,-id_col]
    id_col <- which(colnames(data_quant) == row_names)
    colnames(data_quant)[id_col] <- 'identifier'

    data_quant <- data_quant |>
      dplyr::filter(!is.na(identifier))
  }

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
    dplyr::group_by(sample, replicate, identifier) %>%
    dplyr::summarise(abundance = sum(abundance),
                     .groups = 'drop') %>%
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
                             show_rownames = row_names_bool,
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
                             main = glue::glue("Heatmap: {row_names} {data$quantitative_source} \n")
  )

  if(!'null device' %in% names(grDevices::dev.cur())) { invisible(grDevices::dev.off()) }

  return(plot_save(plot,
                   data,
                   glue::glue("{data$analyte}_{data$quantitative_source}_heatmap"),
                   ...))
}
