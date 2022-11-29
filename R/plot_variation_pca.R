#' Plot the PCA variation in normalized values
#'
#' @description
#' `plot_variation_pca()` is a GGplot2 implementation for plotting the variability in
#' normalized values by PCA analysis, generating two facets. The left facet is a plot of CVs for
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
#' ecoli_proteins %>%
#'   normalize(.method = c("linear", "loess", "randomforest")) %>%
#'   plot_variation_pca()
#'
plot_variation_pca <- function(
    data = NULL,
    ...
){

  # visible bindings
  origin <- NULL
  identifier <- NULL
  abundance_log2 <- NULL
  princom <- NULL
  values <- NULL
  var <- NULL
  val <- NULL
  pcn <- NULL

  check_data(data)

  analyte <- data$analyte

  quant_values <- data %>% get_quant_names()
  data_quant <- data %>% extract(values = quant_values)

  if(!c('origin') %in% names(data_quant)) {
    cli::cli_abort(c("x" = "Input {.emph data} has not been normalized to compare methods"))
  }

  tb_quant <- data_quant %>% transform_log2() %>%
    tidyr::unite(sample, sample, replicate)

  pca_sum <- list()
  for( value in quant_values ){
    df_pca <- tb_quant %>%
      dplyr::filter(origin == value) %>%
      dplyr::select(identifier, sample, abundance_log2) %>%
      tidyr::pivot_wider(
        names_from = "identifier",
        values_from = "abundance_log2",
        values_fill = 0
      ) %>%
      tibble::column_to_rownames('sample') %>%
      as.data.frame()

    df_pca <- df_pca[,which(apply(df_pca, 2, var, na.rm=TRUE) != 0)]

    pca_res <- stats::prcomp(df_pca, scale. = TRUE)
    pca_sum[[value]] <- base::summary(pca_res)$importance %>%
      tidyr::as_tibble() %>%
      dplyr::mutate(var = c('stdev','prop_var','cum_var')) %>%
      tidyr::pivot_longer(!var, names_to = 'princom', values_to = 'val') %>%
      dplyr::mutate(pcn = sub("[A-Z]+", "", princom) %>% as.numeric(),
                    values = value)
  }

  pca_sum <- pca_sum %>%
    dplyr::bind_rows() %>%
    dplyr::filter(var == 'prop_var',
                  val != 0) %>%
    dplyr::group_by(values) %>%
    dplyr::mutate(val = cumsum(val)) %>%
    dplyr::ungroup()

  plot <- pca_sum %>%
    dplyr::mutate(values = forcats::fct_relevel(values,
                                                intersect(
                                                  c('raw', 'median', 'linear', 'limma', 'loess', 'randomforest'),
                                                  unique(pca_sum$values)))) %>%
    ggplot2::ggplot(ggplot2::aes(pcn, val)) +
    ggplot2::geom_point(ggplot2::aes(color=values), alpha=0.5) +
    ggplot2::geom_line(ggplot2::aes(color=values)) +
    ggplot2::scale_x_continuous(breaks = unique(pca_sum$pcn)) +
    ggplot2::scale_y_continuous(limits = c(0,ceiling(max(pca_sum$val)*10)/10)) +
    ggplot2::labs(x = "Principal Component",
                  y = "Cumulative Variance",
                  title ="PCA Variance ~ Normalization",
                  subtitle = "protein abundance log2 co-variance") +
    ggplot2::scale_color_manual(values = theme_palette()) +
    ggplot2::theme_minimal()

  return(plot_save(plot,
                   data,
                   glue::glue("{data$analyte}_normalized_variation_pca"),
                   ...))
}

