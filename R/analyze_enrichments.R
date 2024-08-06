#' Analysis tables and plots of expression values
#'
#' @description
#' `analyze_enrichments()` is a GGplot2 implementation for plotting the expression differences
#' as foldchange ~ statistical significance. See also `plot_proportion()`.  This function can
#' take either a tidyproteomics data object or a table with the required headers.
#'
#' @param data a character defining the column name of the log2 foldchange values.
#' @param top_n a numerical value defining the number of terms to display in the plot
#' @param significance_max a numeric defining the maximum statistical significance to highlight.
#' @param enriched_up_color a color to assign the up enriched values
#' @param enriched_down_color a color to assign the down enriched values
#' @param width a numeric

#'
#' @return a tidyproteomics data object
#' @export
#'
#' @examples
#' library(dplyr, warn.conflicts = FALSE)
#' library(tidyproteomics)
#' hela_proteins %>%
#'    expression(knockdown/control) %>%
#'    analyze_expressions(log2fc_min = 0.5, significance_column = "p_value")
#'
analyze_enrichments <- function(
    data = NULL,
    top_n = 50,
    significance_max = 0.05,
    enriched_up_color = 'blue',
    enriched_down_color = 'red',
    height = 6.5,
    width = 10
){

  # visible bindings
  metric <- NULL

  list.append <- function (x, i) {
    x[[length(x) + 1]] <- i
    x
  }

  if(!is.numeric(significance_max) | significance_max < 0 | significance_max > 1){
    cli::cli_abort("significance_max must be a numeric between 0-1")
  }
  if(!is.numeric(top_n) | top_n < 0 | top_n > 100){
    cli::cli_abort("significance_max must be a numeric between 0-100")
  }

  check_data(data)

  for( set_expression in names(data$analysis)){

    if(is.null(data$analysis[[set_expression]]$enrichment)){ next }
    tbl_enrich <- list()

    for( term in names(data$analysis[[set_expression]]$enrichment) ){
      tbl_enrich <- tbl_enrich %>%
        list.append(
          data$analysis[[set_expression]]$enrichment[[term]]$data %>%
            dplyr::mutate(term = term) %>%
            dplyr::relocate(term)
        )
    }

    tbl_enrich <- tbl_enrich %>% dplyr::bind_rows()
    n_terms <- tbl_enrich %>% nrow()
    plt_title <- set_expression %>%
      stringr::str_replace("/", " / ") %>%
      stringr::str_replace_all("_", " ") %>%
      stringr::str_to_title()

    plt <- tbl_enrich %>%
      dplyr::slice_min(adj_p_value, n = top_n, with_ties = FALSE) %>%
      dplyr::mutate(is_significant = adj_p_value <= significance_max,
                    is_enriched = enrichment > 0) %>%
      dplyr::mutate(term = term %>% stringr::str_replace_all("gene_ontology\\_", "gene_ontology\n")) %>%
      dplyr::mutate(annotation = forcats::fct_reorder(annotation, -enrichment)) %>%
      ggplot2::ggplot(ggplot2::aes(enrichment, annotation)) +
      # ggplot2::geom_col(ggplot2::aes(color = is_enriched), fill = NA) +
      ggplot2::geom_col(ggplot2::aes(fill = is_enriched, alpha = is_significant)) +
      ggplot2::facet_grid(term~., space = 'free', scales = 'free', switch = 'y') +
      ggplot2::scale_fill_manual(values = c('TRUE' = enriched_up_color, 'FALSE' = enriched_down_color)) +
      ggplot2::scale_color_brewer(palette = 'Set1') +
      ggplot2::scale_alpha_manual(values = c(.25, 1)) +
      ggplot2::theme_bw() +
      ggplot2::scale_y_discrete(position = "right") +
      ggplot2::theme(strip.text.y.left = ggplot2::element_text(angle = 0),
                     legend.position = 'none') +
      ggplot2::xlab("Normalized Enrichment Value") +
      ggplot2::ylab("") +
      ggplot2::labs(title = plt_title,
                    subtitle = glue::glue("{n_terms} terms tested; top {top_n} shown; bold = adj_p_value <= {significance_max}"))

    # plot the expressions
    plt_name <- glue::glue("{data$analyte}_enrichment_{stringr::str_replace_all(set_expression, '/', '-')}.png")
    ggplot2::ggsave(filename = plt_name, plot = plt, height = height, width = width)

    # save the expression table
    data_name <- glue::glue("table_{data$analyte}_enrichment_{stringr::str_replace_all(set_expression, '/', '-')}.csv")
    tbl_enrich %>% readr::write_csv(data_name)

  }

  return(data)
}
