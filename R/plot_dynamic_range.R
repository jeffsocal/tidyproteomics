#' Plot CVs by abundance
#'
#' @description
#' `plot_dynamic_range()` is a GGplot2 implementation for plotting the normalization
#' effects on CVs by abundance, visualized as a 2d density plot. Layered on top
#' is a loess smoothed regression of the CVs by abundance, with the median CV
#' shown in _red_ and the dynamic range represented as a box plot on top. The
#' point of this plot is to examine how CVs were minimized through out the abundance
#' profile. Some normalization methods function well at high abundance yet leave
#' retain high CVs at lower abundance.
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
#' hela_proteins %>%
#'   normalize(.method = c("linear", "loess", "randomforest")) %>%
#'   plot_dynamic_range()
#'
plot_dynamic_range <- function(
    data = NULL,
    ...
){

  # visible bindings
  abundance <- NULL
  origin <- NULL
  identifier <- NULL
  abundance_sd <- NULL
  abundance_mean <- NULL
  abundance_cv <- NULL
  abundance_min <- NULL
  abundance_max <- NULL
  abundance_95l <- NULL
  abundance_95h <- NULL
  cv_mean <- NULL

  check_data(data)

  quant_values_all <- c('raw', 'median', 'scaled', 'linear',
                        'limma', 'loess', 'svm', 'randomforest')
  quant_values <- get_quant_names(data)
  data_quant <- data %>% extract(values = quant_values)

  data_quant_norm_range <- data_quant %>%
    dplyr::mutate(abundance = abundance %>% log10(),
                  origin = origin %>% as.factor()) %>%
    dplyr::filter(!is.infinite(abundance)) %>%
    dplyr::group_by(origin, sample) %>%
    dplyr::summarise(
      abundance_min = min(abundance, na.rm = T),
      abundance_max = max(abundance, na.rm = T),
      abundance_95l = stats::quantile(abundance, .025, na.rm=T),
      abundance_95h = stats::quantile(abundance, .975, na.rm=T),
      .groups = "drop"
    ) %>%
    dplyr::ungroup()

  data_quant_cvs <- data_quant %>%
    dplyr::mutate(origin = sub("data_vals_", "", origin),
                  origin = sub("\\.rds", "", origin),
                  origin = origin %>% as.factor()) %>%
    dplyr::group_by(origin, identifier, sample) %>%
    dplyr::summarise(
      abundance_mean = mean(abundance, na.rm = T),
      abundance_sd = stats::sd(abundance, na.rm = T),
      abundance_cv = abundance_sd / abundance_mean,
      .groups = "drop"
    ) %>%
    dplyr::ungroup() %>%
    dplyr::filter(!is.na(abundance_cv))

  data_quant_cvs_mean <- data_quant_cvs %>%
    dplyr::group_by(origin, sample) %>%
    dplyr::summarise(
      cv_mean = mean(abundance_cv, na.rm = T),
      cv_sd = stats::sd(abundance_cv, na.rm = T),
      .groups = "drop"
    )


  range_y <- max(data_quant_cvs$abundance_cv) * 1.1
  range_x <- stats::median(data_quant_cvs$abundance_mean %>% log10())
  min_x <- min(data_quant_cvs$abundance_mean %>% log10())

  plot <- data_quant_cvs %>%
    dplyr::mutate(origin = forcats::fct_relevel(origin,
                                                intersect(
                                                  quant_values_all,
                                                  unique(data_quant_cvs$origin)))) %>%
    ggplot2::ggplot(ggplot2::aes(log10(abundance_mean), abundance_cv)) +
    ggplot2::geom_hex() +
    ggplot2::geom_smooth(fill=NA, color='red', method='gam',
                         formula = stats::as.formula('y ~ s(x, bs = "cs")')) +
    ggplot2::geom_segment(data=data_quant_norm_range,
                          ggplot2::aes(x=abundance_min, y=range_y,
                                       xend=abundance_max, yend=range_y),
                          color='dodgerblue') +
    ggplot2::geom_segment(data=data_quant_norm_range,
                          ggplot2::aes(x=abundance_95l, y=range_y,
                                       xend=abundance_95h, yend=range_y),
                          color='dodgerblue', size=2) +
    ggplot2::geom_hline(yintercept = range_y*1.2, color=NA) +
    ggplot2::geom_point(ggplot2::aes(x=range_x, y=range_y), color='lightblue') +
    ggplot2::geom_text(data=data_quant_norm_range,
                       ggplot2::aes(x=abundance_95l, y=range_y,
                                    label=paste('range', signif(abundance_max-abundance_min, 2))),
                       size=2.5, vjust=-0.7, hjust=0) +
    ggplot2::geom_text(data=data_quant_cvs_mean,
                       ggplot2::aes(x=min_x, y=range_y,
                                    label=signif(cv_mean, 2)),
                       size=3.5, vjust=2, hjust=0,
                       color='red') +
    ggplot2::scale_fill_viridis_c() +
    ggplot2::facet_grid(sample~origin) +
    ggplot2::scale_color_manual(values = theme_palette()) +
    ggplot2::theme_bw() +
    ggplot2::labs(
      title = "Normalization Effect on Variation",
      subtitle = "CVs (red) ~ Abundance (dynamic range)",
      x = 'Abundance log10',
      y = 'Abundance CV')

  return(plot_save(plot,
                   data,
                   glue::glue("{data$analyte}_dynamic_range"),
                   ...))
}
