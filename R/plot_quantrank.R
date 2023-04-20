#' Plot the variation in normalized values
#'
#' @description
#' `plot_quantrank()` is a GGplot2 implementation for plotting the variability in
#' normalized values, generating two facets. The left facet is a plot of CVs for
#' each normalization method. The right facet is a plot of the 95%CI in abundance,
#' essentially the conservative dynamic range. The goal is to select a normalization
#' method that minimizes CVs while also retaining the dynamic range.
#'
#' @param data tidyproteomics data object
#' @param accounting character string
#' @param type character string
#' @param show_error a boolean
#' @param show_rank_scale a boolean
#' @param limit_rank a numerical vector of 2
#' @param display_filter a numeric between 0 and 1
#' @param display_cutoff a numeric between 0 and 1
#' @param palette a string representing the palette for scale_fill_brewer()
#' @param impute_max a numeric representing the largest allowable imputation percentage
#' @param ... passthrough for ggsave see `plotting`
#'
#' @return a (tidyproteomics data-object | ggplot-object)
#' @export
#'
#' @examples
#' library(dplyr, warn.conflicts = FALSE)
#' library(tidyproteomics)
#' hela_proteins %>% plot_quantrank()
#'
#' hela_proteins %>% plot_quantrank(type = "lines")
#'
#' hela_proteins %>% plot_quantrank(display_filter = "log2_foldchange", display_cutoff = 1)
#'
#' hela_proteins %>% plot_quantrank(limit_rank = c(1,50), show_rank_scale = TRUE)
#'
plot_quantrank <- function(
    data = NULL,
    accounting = 'protein',
    type = c("points", "lines"),
    show_error = TRUE,
    show_rank_scale = FALSE,
    limit_rank = NULL,
    display_subset = NULL,
    display_filter = c('none','log2_foldchange','p_value','adj_p_value'),
    display_cutoff = 1,
    palette = 'YlGnBu',
    impute_max = 0.5,
    ...
){

  # visible bindings
  log2_foldchange <- NULL
  p_value <- NULL
  identifier <- NULL
  adj_p_value <- NULL
  abundance <- NULL
  abundance_median <- NULL
  disp_filter <- NULL
  abundance_975 <- NULL
  abundance_025 <- NULL

  check_data(data)
  v_terms <- data %>% get_annotation_terms()

  type <- rlang::arg_match(type)
  display_filter <- rlang::arg_match(display_filter)
  accounting <- rlang::arg_match(accounting, v_terms)
  if(!is.logical(show_error)) {cli::cli_abort("show_error needs to be TRUE or FALSE, not `{show_error}`")}
  if(!is.logical(show_rank_scale)) {cli::cli_abort("show_rank_scale needs to be TRUE or FALSE, not `{show_rank_scale}`")}
  if(!is.null(limit_rank)) {
    if(!is.numeric(limit_rank) || length(limit_rank) != 2) {cli::cli_abort("limit_rank needs to be vector of 2 numerics")}
  }

  tb_summary <- data %>% table_quantrank(accounting, display_filter)
  tb_med <- tb_summary$median
  tb_sum <- tb_sum_org <- tb_summary$summary

  r_min <- min(tb_med$rank)
  r_max <- max(tb_med$rank)
  if(!is.null(limit_rank)) {
    r_min <- min(limit_rank)
    r_max <- max(limit_rank)
  }

  tb_sum <- tb_sum_org <- tb_sum %>% dplyr::filter(rank >= r_min) %>% dplyr::filter(rank <= r_max)
  tb_med <- tb_med %>% dplyr::filter(rank >= r_min) %>% dplyr::filter(rank <= r_max)

  # y_min <- 10 ^ floor(log10(min(tb_med$abundance_median)))
  # y_max <- 10 ^ ceiling(log10(max(tb_med$abundance_median)))

  if(display_filter != 'none') {type <- 'points'}

  if(type == 'points'){
    plot <- tb_med %>%
      dplyr::rename(abundance = abundance_median) %>%
      ggplot2::ggplot(ggplot2::aes(identifier, abundance))
  } else {
    plot <- tb_med %>%
      dplyr::rename(abundance = abundance_median) %>%
      ggplot2::ggplot(ggplot2::aes(rank, abundance))
  }

  for( i in 1:2 ){

    tb_sum <- tb_sum_org
    if(display_filter == 'none' & is.null(display_subset)){ i <- 2 }

    if(i == 1) {
      tb_sum <- tb_sum %>% dplyr::mutate(sample = NA)
      alpha = .1
    }
    if(i == 2) {
      if(display_filter != 'none'){
        tb_sum$disp_filter <- unlist(tb_sum[,display_filter])
        if(grepl("p_value", display_filter)) {
          tb_sum <- tb_sum %>% dplyr::filter(disp_filter <= display_cutoff)
        } else {
          tb_sum <- tb_sum %>% dplyr::filter(disp_filter > display_cutoff)
        }
      } else if(!is.null(display_subset)){
        dsp <- display_subset %>% extract() %>% dplyr::select(identifier) %>% unique()
        tb_sum <- tb_sum %>% dplyr::inner_join(dsp, by = 'identifier')
        display_filter <- "subset"
      }
      alpha = .25
    }

    if(type == 'points'){

      if(show_error == TRUE){
        plot <- plot +
          ggplot2::geom_segment(data = tb_sum,
                                ggplot2::aes(identifier, abundance_975,
                                             xend=identifier, yend=abundance_025,
                                             color=sample), size=0.5, alpha=alpha)
      }
      plot <- plot +
        ggplot2::geom_point(data = tb_sum,
                            ggplot2::aes(identifier, abundance_median,
                                         color=sample), size=1, alpha=alpha*2)
    } else {
      plot <- tb_med %>%
        dplyr::rename(abundance = abundance_median) %>%
        ggplot2::ggplot(ggplot2::aes(rank, abundance_median)) +
        ggplot2::geom_smooth(data = tb_sum, formula = stats::as.formula('y ~ s(x, bs = "cs")'),
                             method = 'gam', fill = NA,
                             ggplot2::aes(rank, abundance_median,
                                          color=sample), linetype = 1, size=alpha*2)

      if(show_error == TRUE){
        plot <- plot +
          ggplot2::geom_smooth(data = tb_sum, formula = stats::as.formula('y ~ s(x, bs = "cs")'),
                               method = 'gam', fill = NA,
                               ggplot2::aes(rank, abundance_975,
                                            color=sample), linetype = 1, size=alpha) +
          ggplot2::geom_smooth(data = tb_sum, formula = stats::as.formula('y ~ s(x, bs = "cs")'),
                               method = 'gam', fill = NA,
                               ggplot2::aes(rank, abundance_025,
                                            color=sample), linetype = 1, size=alpha)
      }
    }
  }

  if(type == 'points' & display_filter != 'none' & i == 2){
    plot <- plot +
      ggrepel::geom_text_repel(data = tb_med %>%
                                 dplyr::rename(abundance = abundance_median) %>%
                                 dplyr::filter(identifier %in% tb_sum$identifier),
                               ggplot2::aes(label = identifier),
                               size = 3, box.padding = 1)
  }

  subtitle <- glue::glue("Normalization: {data$quantitative_source}")
  if(display_filter != 'none') {
    subtitle <- glue::glue("{subtitle} Filter: {display_filter} {display_cutoff}")
  }

  plot <- plot +
    ggplot2::scale_y_log10() +
    ggplot2::theme_classic() +
    ggplot2::labs(title = glue::glue("Quantitation Rank Plot: {data$origin}"),
                  subtitle = subtitle,
                  x = glue::glue('{accounting} rank'), y = glue::glue("log10 abundance")) +
    ggplot2::scale_color_manual(values = theme_palette()) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90))

  if(show_rank_scale == FALSE){
    plot <- plot +
      ggplot2::theme(axis.ticks.x = ggplot2::element_blank(),
                     axis.text.x = ggplot2::element_blank(),
                     panel.grid.major.x = ggplot2::element_blank())
  }

  return(plot_save(plot,
                   data,
                   glue::glue("{data$analyte}_quantitation_rank"),
                   ...))
}



#' Helper function to quantitation plots
#'
#' @description
#' `table_quantrank()`
#'
#' @param data tidyproteomics data object
#' @param accounting character string
#' @param display_filter a numeric between 0 and 1
#'
#' @return a (tidyproteomics data-object | ggplot-object)
#'
#' @examples
#' library(dplyr, warn.conflicts = FALSE)
#' library(tidyproteomics)
#' hela_proteins %>% plot_quantrank()
#'
#' hela_proteins %>% plot_quantrank(type = 'lines')
#'
#' hela_proteins %>% plot_quantrank(type = 'lines', display_filter = 'log2_foldchange', display_cutoff = 1)
#'
table_quantrank <- function(
    data = NULL,
    accounting = 'protein',
    display_filter = c('none','log2_foldchange','p_value','adj_p_value')
){

  # visible bindings
  log2_foldchange <- NULL
  p_value <- NULL
  identifier <- NULL
  adj_p_value <- NULL
  abundance <- NULL
  abundance_median <- NULL
  disp_filter <- NULL
  abundance_975 <- NULL
  abundance_025 <- NULL

  check_data(data)
  v_terms <- data %>% get_annotation_terms()

  display_filter <- rlang::arg_match(display_filter)
  accounting <- rlang::arg_match(accounting, v_terms)
  # if(deviation_filter < 0 || deviation_filter > 1) {cli::cli_abort("x" = "deviation_filter is between 0 and 1, not `{deviation_filter}`")}

  if(display_filter != 'none'){
    pairs <- data %>% get_sample_names() %>% utils::combn(2)
    tb_exp <- list()
    for(i in 1:ncol(pairs)){
      tb_exp[[i]] <- data %>% expression_test(experiment = pairs[1,i], control = pairs[2,i])
    }
    tb_exp <- tb_exp %>%
      dplyr::bind_rows() %>%
      munge_identifier(munge = 'combine', data$identifier) %>%
      dplyr::mutate(log2_foldchange = abs(log2_foldchange)) %>%
      dplyr::mutate(adj_p_value = stats::p.adjust(p_value)) %>%
      dplyr::group_by(identifier) %>%
      dplyr::summarise(
        log2_foldchange = max(log2_foldchange, na.rm = TRUE),
        p_value = min(p_value, na.rm = TRUE),
        adj_p_value = min(adj_p_value, na.rm = TRUE),
        .groups = 'drop'
      )
  }

  tb <- data %>% meld(single_quant_source = TRUE) %>%
    dplyr::filter(!is.na(abundance)) %>%
    dplyr::filter(abundance > 0) %>%
    dplyr::select(dplyr::matches(glue::glue("^{accounting}$|abundance|sample|replicate")))

  tb$identifier <- unlist(tb[,accounting])

  tb_med <- tb %>%
    dplyr::group_by(identifier) %>%
    dplyr::summarise(abundance_median = stats::median(abundance, na.rm=TRUE),
                     .groups = 'drop') %>%
    dplyr::arrange(dplyr::desc(abundance_median)) %>%
    dplyr::mutate(rank = dplyr::row_number()) %>%
    dplyr::mutate(identifier = forcats::fct_reorder(identifier, rank))

  if(display_filter != 'none') {tb_med <- tb_med %>% dplyr::left_join(tb_exp, by = 'identifier')}

  tb_sum <- tb %>%
    dplyr::group_by(identifier, sample) %>%
    dplyr::summarise(abundance_median = stats::median(abundance, na.rm=TRUE),
                     abundance_975 = stats::quantile(abundance, 0.975, na.rm=TRUE),
                     abundance_025 = stats::quantile(abundance, 0.025, na.rm=TRUE),
                     .groups = 'drop') %>%
    dplyr::inner_join(tb_med %>% dplyr::select(!dplyr::matches('median')), by = c("identifier")) %>%
    dplyr::mutate(identifier = forcats::fct_reorder(identifier, rank))

  return(list(median = tb_med, summary = tb_sum))

}
