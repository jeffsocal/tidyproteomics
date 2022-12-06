#' Plot PCA values
#'
#' @description
#' `plot_pca()` is a GGplot2 implementation for plotting two principal components
#' from a PCA analysis, visualized as a scatter.
#'
#' @param data tidyproteomics data object
#' @param variables a character vector of the 2 PCs to plot. Acceptable values
#' include (PC1, PC2, PC3 ... PC9). Default c('PC1','PC2').
#' @param labels a boolean
#' @param label_size a numeric
#' @param ... passthrough for ggsave see `plotting`
#'
#' @return a (tidyproteomics data-object | ggplot-object)
#' @export
#'
#' @examples
#' library(dplyr, warn.conflicts = FALSE)
#' library(tidyproteomics)
#' hela_proteins <- hela_proteins %>%
#'   normalize(.method = c("scaled", "median", "linear", "limma", "loess")) %>%
#'   select_normalization()
#'
#' hela_proteins %>% plot_pca()
#'
#' # a different PC set
#' hela_proteins %>% plot_pca(variables = c("PC2", "PC3"))
#'
#' # a PC scree plot
#' hela_proteins %>% plot_pca("scree")
#'
plot_pca <- function(
    data = NULL,
    variables = c('PC1','PC2'),
    labels = TRUE,
    label_size = 3,
    ...
){

  # visible bindings
  identifier <- NULL
  abundance <- NULL
  var <- NULL
  princom <- NULL
  val <- NULL
  pcn <- NULL
  sample_exp <- NULL

  variables <- variables[1:min(length(variables), 2)]
  variables <- rlang::arg_match(variables, c('scree', paste0("PC", 1:9)), multiple = T)
  check_data(data)

  analyte <- data$analyte
  quantval <- data$quantitative_source
  data_quant <- data %>% extract(values = quantval, na.rm = TRUE)

  dt_pca <- data_quant %>%
    dplyr::mutate(sample = paste(sample,replicate, sep=".")) %>%
    dplyr::select(identifier, sample, abundance) %>%
    tidyr::pivot_wider(
      names_from = 'identifier',
      values_from = 'abundance',
      values_fill = 0
    ) %>%
    tibble::column_to_rownames('sample') %>%
    as.data.frame()

  dt_pca <- dt_pca[,which(apply(dt_pca, 2, var, na.rm=TRUE) != 0)]

  pca_res <- stats::prcomp(dt_pca, scale. = TRUE)
  pca_sum <- base::summary(pca_res)$importance %>%
    tibble::as_tibble() %>%
    dplyr::mutate(var = c('stdev','prop_var','cum_var')) %>%
    tidyr::pivot_longer(!var, names_to = 'princom', values_to = 'val') %>%
    dplyr::mutate(pcn = sub("[A-Z]+", "", princom) %>% as.numeric())

  if(length(variables) == 1 && variables == 'scree'){

    plot <- pca_sum %>%
      dplyr::mutate(val = val * 100) %>%
      dplyr::filter(var == 'prop_var') %>%
      ggplot2::ggplot(ggplot2::aes(pcn, val)) +
      ggplot2::geom_bar(stat='identity') +
      ggplot2::geom_text(ggplot2::aes(label = paste0(signif(val, 2), "%")), hjust=0.5, vjust=-0.2) +
      ggplot2::scale_x_continuous(breaks = unique(pca_sum$pcn)) +
      ggplot2::scale_y_continuous(limits = c(0,100)) +
      ggplot2::labs(x = "Principal Component",
                    y = "Proportion of Total Variance",
                    title = "PCA Scree Plot",
                    subtitle = "protein abundance log2 co-variance") +
      ggplot2::theme_classic()

  } else {

    pcas <- pca_sum %>%
      dplyr::select(princom) %>%
      unique() %>%
      unlist() %>%
      as.character()

    pcas <- intersect(pcas, variables[1:2])

    pcx <- pcas[1]
    pcy <- pcas[2]
    pcx_p <- pca_sum %>%
      dplyr::filter(princom == pcx & var == 'prop_var') %>%
      dplyr::select(val) %>% unlist() %>% as.numeric()

    pcy_p <- pca_sum %>%
      dplyr::filter(princom == pcy & var == 'prop_var') %>%
      dplyr::select(val) %>% unlist() %>% as.numeric()

    plot <- pca_res$x %>%
      as.data.frame() %>%
      dplyr::rename(pcx = dplyr::all_of(pcx),
                    pcy = dplyr::all_of(pcy)) %>%
      tibble::rownames_to_column('sample_exp') %>%
      tidyr::separate(sample_exp, into = c('sample','replicate'), sep="\\.", remove = F) %>%
      ggplot2::ggplot(ggplot2::aes(pcx, pcy)) +
      ggplot2::geom_point(ggplot2::aes(color=sample), alpha=0.5, size=5)

    if(labels == TRUE) {
      plot <- plot +
        ggrepel::geom_text_repel(ggplot2::aes(label = sample_exp), size = label_size)
    }

    plot <- plot +
      ggplot2::scale_color_manual(values = theme_palette()) +
      ggplot2::labs(x = paste0(pcx, " (",signif(pcx_p * 100, 3),"%)"),
                    y = paste0(pcy, " (",signif(pcy_p * 100, 3),"%)"),
                    title = paste("PCA analysis", pcx, "~", pcy),
                    subtitle = paste0("protein abundance log2 co-variance (", quantval, ")")) +
      ggplot2::theme_classic() +
      ggplot2::theme(axis.text = ggplot2::element_blank(),
                     axis.ticks = ggplot2::element_blank())
  }

  return(plot_save(plot,
                   data,
                   glue::glue("{data$analyte}_{data$quantitative_source}_pca_{paste(variables, collapse='')}"),
                   ...))
}

