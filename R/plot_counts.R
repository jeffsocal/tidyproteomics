#' Plot the variation in normalized values
#'
#' @description
#' `plot_counts()` is a GGplot2 implementation for plotting counting statistics.
#'
#' @param data tidyproteomics data object
#' @param accounting character string
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
plot_counts <- function(
    data = NULL,
    accounting = c('protein_groups','proteins','peptides','peptides_unique'),
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
  match_between_runs <- NULL
  imputed <- NULL

  accounting <- rlang::arg_match(accounting)
  check_data(data)
  limit <- nrow(data$experiments)

  fat <- list(
    data %>%
      subset(match_between_runs == FALSE, .verbose = FALSE) %>%
      subset(imputed >= !!impute_max, .verbose = FALSE) %>%
      summary('sample_id', destination = 'return', limit = limit) %>% dplyr::mutate(mbr = 'xMBR'),
    data %>% summary('sample_id', destination = 'return', limit = limit) %>% dplyr::mutate(mbr = 'MBR')
  )  %>%
    dplyr::bind_rows() %>%
    dplyr::full_join(
      data$experiment %>%
        dplyr::select(c('sample_id', 'sample', 'replicate')),
      by = 'sample_id'
    )

  w <- which(duplicated(fat[,c('sample_id', 'sample', accounting)]))
  if(length(w) > 0) {fat <- fat[-w,]}

  fat$metric <- unlist(fat[,accounting])
  legend_title <- 'Replicate'
  fat$fill <- fat$replicate
  if(length(unique(fat$replicate)) > 9) {
    fat$fill <- fat$sample
    legend_title <- 'Sample'
  }
  metric_max <- max(fat$metric) * 1.1

  fat_mean <- fat %>%
    dplyr::group_by(sample, mbr) %>%
    dplyr::summarise(
      replicate = mean(as.numeric(replicate)),
      q025 = round(stats::quantile(metric, .025)),
      q975 = round(stats::quantile(metric, .975)),
      ci95 = q975 - q025,
      metric = round(stats::median(metric)),
      .groups = 'drop') %>%
    dplyr::mutate(label = ifelse(ci95 == 0, metric,
                                 glue::glue("{metric} \u00B1 {ci95}")),
                  label = ifelse(mbr == 'MBR', glue::glue("{label} mbr"), label),
                  metric = ifelse(mbr == 'MBR', metric * 1.05, metric * .95))

  plot <- fat %>%
    # dplyr::mutate(replicate = as.factor(replicate)) %>%
    ggplot2::ggplot(ggplot2::aes(replicate, metric)) +
    ggplot2::geom_bar(ggplot2::aes(group = replicate), color = 'grey',
                      stat = 'identity', fill=NA, position = 'identity',
                      width=1) +
    ggplot2::geom_bar(ggplot2::aes(fill=fill, group = replicate),
                      stat = 'identity', alpha=.67, position = 'identity',
                      width=1) +
    ggplot2::facet_grid(.~sample, scales = 'free_x', space = 'free_x')  +
    ggplot2::geom_text(data = fat_mean,
                       ggplot2::aes(label = label),
                       size = 3) +
    ggplot2::geom_hline(yintercept = metric_max, color=NA) +
    ggplot2::scale_fill_brewer(palette = palette) +
    ggplot2::scale_color_brewer(palette = palette) +
    ggplot2::theme_classic() +
    ggplot2::labs(title = paste(data$origin, "--", sub("s$", "", accounting), 'counts'),
                  fill = legend_title) +
    ggplot2::ylab(accounting) +
    ggplot2::theme(axis.title.x = ggplot2::element_blank(),
                   axis.ticks.x = ggplot2::element_blank(),
                   axis.text.x = ggplot2::element_blank(),
                   panel.grid.major.x = ggplot2::element_blank())

  return(plot_save(plot,
                   data,
                   glue::glue("{data$analyte}_counts"),
                   ...))
}
