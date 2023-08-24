#' Plot the variation in normalized values
#'
#' @description
#' `plot_variation_cv()` is a GGplot2 implementation for plotting the variability in
#' normalized values, generating two facets. The left facet is a plot of CVs for
#' each normalization method. The right facet is a plot of the 95%CI in abundance,
#' essentially the conservative dynamic range. The goal is to select a normalization
#' method that minimizes CVs while also retaining the dynamic range.
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
#'   normalize(.method = c("scaled", "median", "linear", "limma", "loess")) %>%
#'   plot_variation_cv()
#'
plot_variation_cv <- function(
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
  cv_mean <- NULL
  abundance_95h <- NULL
  abundance_95l <- NULL
  dynamic_range <- NULL

  check_data(data)

  analyte <- data$analyte

  quant_values <- get_quant_names(data)
  data_quant <- data %>% extract(values = quant_values)

  quant_values_all <- c('raw', 'median', 'scaled', 'linear',
                        'limma', 'loess', 'svm', 'randomforest')

  norm_vals <- intersect(quant_values_all, quant_values)
  if(length(norm_vals) == 1) {
    cli::cli_alert_warning("Normalization has not yet been performed.")
    return(data)
  }

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
      dynamic_range = abundance_95h - abundance_95l,
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

  n_colors <- data_quant_cvs_mean$sample %>% unique() %>% length()

  p_cv <- data_quant_cvs_mean %>%
    dplyr::mutate(origin = forcats::fct_relevel(origin, norm_vals)) %>%
    ggplot2::ggplot(ggplot2::aes(origin, cv_mean, color=sample)) +
    ggplot2::geom_point() +
    ggplot2::geom_line(ggplot2::aes(group=sample)) +
    ggplot2::geom_hline(yintercept = 0, color=NA) +
    ggplot2::scale_y_continuous(n.breaks = 7) +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle=20, hjust=1),
                   legend.position = 'none',
                   axis.title.y = ggplot2::element_text(size=9)) +
    ggplot2::scale_color_manual(values = theme_palette(n_colors)) +
    ggplot2::labs(subtitle = "Quantitative Variation") +
    ggplot2::xlab('') +
    ggplot2::ylab('CVs (sd/mean)')

  p_dr <- data_quant_norm_range %>%
    dplyr::mutate(origin = forcats::fct_relevel(origin, norm_vals)) %>%
    ggplot2::ggplot(ggplot2::aes(origin, dynamic_range, color=sample)) +
    ggplot2::geom_point() +
    ggplot2::geom_line(ggplot2::aes(group=sample)) +
    ggplot2::scale_y_continuous(n.breaks = 11) +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle=20, hjust=1),
                   axis.title.y = ggplot2::element_text(size=9)) +
    ggplot2::scale_color_manual(values = theme_palette(n_colors)) +
    ggplot2::labs(subtitle = "Quantitative Dynamic Range") +
    ggplot2::xlab('') +
    ggplot2::ylab('Dynamic Range (95%CI Log10)')

  # plot <- data_quant_cvs_mean %>%
  #   dplyr::select(origin, sample, mean = cv_mean) %>%
  #   dplyr::mutate(stat = 'CVs (sd/mean)') %>%
  #   dplyr::bind_rows(
  #     data_quant_norm_range %>%
  #       dplyr::mutate(mean = abundance_95h - abundance_95l) %>%
  #       dplyr::select(origin, sample, mean) %>%
  #       dplyr::mutate(stat = 'Dynamic Range (95%CI Log10)')
  #   ) %>%
  #   dplyr::mutate(origin = forcats::fct_relevel(origin, norm_vals)) %>%
  #   ggplot2::ggplot(ggplot2::aes(origin, mean, color=sample)) +
  #   ggplot2::geom_point() +
  #   ggplot2::geom_line(ggplot2::aes(group=sample)) +
  #   ggplot2::geom_hline(yintercept = 0, color=NA) +
  #   ggplot2::scale_y_continuous(n.breaks = 11) +
  #   ggplot2::theme_minimal() +
  #   ggplot2::facet_wrap(~stat, scales='free', strip.position = 'left') +
  #   ggplot2::theme(axis.text.x = ggplot2::element_text(angle=30, hjust=1)) +
  #   ggplot2::scale_color_manual(values = theme_palette()) +
  #   ggplot2::labs(title = "Normalization",
  #                 subtitle = "Effect on Variation and Range",
  #                 x = '', y ='')

  plot <- gridExtra::grid.arrange(p_cv, p_dr, nrow = 1,
                                  widths = c(.85,1.15),
                                  top=grid::textGrob("Normalization Effects on Variation and Range",
                                                     x = 0.05, hjust = 0))

  return(plot_save(plot,
                   data,
                   glue::glue("{data$analyte}_normalizated_variation"),
                   ...))
}
