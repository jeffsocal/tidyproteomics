#' Plot the accounting of proteins. peptides, and other counts
#'
#' @description
#' `plot_counts()` is a GGplot2 implementation for plotting counting statistics.
#'
#' @param data tidyproteomics data object
#' @param accounting character string
#' @param show_replicates boolean to visualize replicates
#' @param impute_max a numeric representing the largest allowable imputation percentage
#' @param palette a string representing the palette for scale_fill_brewer()
#' @param ... passthrough for ggsave see `plotting`
#'
#' @return a (tidyproteomics data-object | ggplot-object)
#' @export
#'
#' @examples
#' library(dplyr, warn.conflicts = FALSE)
#' library(tidyproteomics)
#' hela_proteins %>% plot_counts()
#'
#' hela_proteins %>% plot_counts(show_replicates = FALSE, palette = 'Blues')
#'
plot_counts <- function(
    data = NULL,
    accounting = NULL,
    show_replicates = TRUE,
    impute_max = 0.5,
    palette = 'YlGnBu',
    ...
){

  # visible bindings
  mbr <- NULL
  metric <- NULL
  q975 <- NULL
  q025 <- NULL
  ci95 <- NULL
  label <- NULL
  imputed <- NULL
  fill <- NULL
  is_mbr <- NULL


  check_data(data)

  if(!is.logical(show_replicates)) { cli::cli_abort("show_replicates must be TRUE|FALSE") }
  # if(!is.logical(facet_groups)) { cli::cli_abort("facet_groups must be TRUE|FALSE") }

  accounting_cols <- c(paste0(data$identifier, c('s', '_groups')), 'peptides','peptides_unique')
  if(is.null(accounting)) { accounting <- accounting_cols[1]}
  accounting <- rlang::arg_match(accounting, accounting_cols)

  fat <- data %>% analysis_counts(impute_max = impute_max)
  n_samples <- data %>% get_sample_names() %>% length()

  w <- which(duplicated(fat[,c('sample_id', 'sample', accounting)]))
  if(length(w) > 0) {fat <- fat[-w,]}

  fat$metric <- unlist(fat[,accounting])
  legend_title <- 'Replicate'
  fat$fill <- fat$replicate
  if(show_replicates == FALSE | length(unique(fat$replicate)) > 9) {
    fat$fill <- fat$sample
    legend_title <- 'Sample'
  } else if('sample_group' %in% colnames(fat)) {
    fat$fill <- fat$sample_group
    legend_title <- 'Group'
  }
  # else if(facet_groups == TRUE & 'sample_group' %in% colnames(fat)) {
  #   fat$fill <- fat$sample_group
  #   legend_title <- 'Group'
  # }

  fat$fill <- fat$fill %>% as.factor()

  metric_max <- max(fat$metric) * 1.1

  fat_mean <- fat %>%
    dplyr::group_by(sample, is_mbr) %>%
    dplyr::summarise(
      replicate = mean(as.numeric(replicate)),
      q025 = round(stats::quantile(metric, .025)),
      q975 = round(stats::quantile(metric, .975)),
      ci95 = q975 - q025,
      metric = round(stats::median(metric)),
      .groups = 'drop') %>%
    dplyr::mutate(label = ifelse(ci95 == 0, metric,
                                 glue::glue("{metric} \u00B1 {ci95}")),
                  label = ifelse(is_mbr == TRUE, glue::glue("{label} mbr"), label),
                  label_pos = ifelse(is_mbr == TRUE, metric * 1.05, metric * .95))


  if(show_replicates == FALSE){ #facet_groups == TRUE){
    plot <- fat_mean %>%
      dplyr::mutate(fill = as.factor(sample)) %>%
      # dplyr::mutate(replicate = as.factor(replicate)) %>%
      ggplot2::ggplot(ggplot2::aes(sample, metric)) +
      ggplot2::theme_minimal() +
      ggplot2::geom_col(color = 'grey', fill=NA,
                        position = 'identity', width=.95) +
      ggplot2::geom_col(ggplot2::aes(fill=fill), alpha=.66,
                        position = 'identity', width=.95) +
      ggplot2::theme(axis.title.x = ggplot2::element_blank(),
                     # axis.ticks.x = ggplot2::element_blank(),
                     # axis.ticks.x = ggplot2::element_line(color='grey'),
                     axis.ticks.y = ggplot2::element_line(color='grey'),
                     legend.position = 'none',
                     panel.grid.major.x = ggplot2::element_blank(),
                     panel.grid.minor.y = ggplot2::element_blank(),
                     panel.grid.major.y = ggplot2::element_blank())

  } else {
    plot <- fat %>%
      # dplyr::mutate(replicate = as.factor(replicate)) %>%
      ggplot2::ggplot(ggplot2::aes(replicate, metric)) +
      ggplot2::theme_classic() +
      ggplot2::geom_col(ggplot2::aes(group = replicate), color = 'grey',
                        fill=NA, position = 'identity', width=1) +
      ggplot2::geom_col(ggplot2::aes(fill=fill, group = replicate),
                        alpha=.66, position = 'identity', width=1) +
      ggplot2::facet_grid(.~sample, scales = 'free_x', space = 'free_x') +
      ggplot2::theme(axis.title.x = ggplot2::element_blank(),
                     axis.ticks.x = ggplot2::element_blank(),
                     axis.text.x = ggplot2::element_blank(),
                     panel.grid.major.x = ggplot2::element_blank())
  }

  plot <- plot +
    # dplyr::mutate(replicate = as.factor(replicate)) %>%
    ggplot2::geom_text(data = fat_mean,
                       ggplot2::aes(y = label_pos, label = label),
                       size = 3) +
    ggplot2::geom_hline(yintercept = metric_max, color=NA) +
    ggplot2::labs(title = paste(data$origin, "--", sub("s$", "", accounting), 'counts'),
                  fill = legend_title) +
    ggplot2::ylab(accounting)

  if(n_samples > 9){
    this_palette <- theme_palette(n_samples)

    plot <- plot +
      ggplot2::scale_fill_manual(values = this_palette) +
      ggplot2::scale_color_manual(values = this_palette)

  } else {

    plot <- plot +
      ggplot2::scale_fill_brewer(palette = palette) +
      ggplot2::scale_color_brewer(palette = palette)
  }


  return(plot_save(plot,
                   data,
                   glue::glue("{data$analyte}_counts"),
                   ...))
}
