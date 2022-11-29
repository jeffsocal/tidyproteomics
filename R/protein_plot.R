#' Visualize mapped sequence data
#'
#' @param mapped_data a tidyproteomics data-object, specifically of sequencing origin
#' @param protein a character string
#' @param row_length a numeric
#' @param samples a character string
#' @param modifications a character string
#' @param ncol a numeric
#' @param nrow a numeric
#' @param color_sequence a character string
#' @param color_modifications a character vector
#' @param show_modification_precent a boolean
#'
#' @return a list of protein mappings
#' @export
#'
#' @examples
#' library(dplyr, warn.conflicts = FALSE)
#' library(tidyproteomics)
#'
#' ecoli_protein_map <- ecoli_peptides %>%
#'    protein_map(fasta = path_to_package_data('fasta'))
#'
#' ecoli_protein_map %>% plot_protein('P0A6Y8')
#'
plot_protein <- function(
    mapped_data = NULL,
    protein = NULL,
    row_length = 50,
    samples = NULL,
    modifications = NULL,
    ncol = NULL,
    nrow = NULL,
    color_sequence = 'grey60',
    color_modifications = c("red","blue", "orange","skyblue","purple","yellow"),
    show_modification_precent = TRUE
){

  # visible bindings
  abundance <- NULL
  modification <- NULL
  frequency <- NULL
  aa <- NULL

  protein_map <- protein_map_munge(mapped_data = mapped_data,
                                   protein = protein,
                                   row_length = row_length,
                                   samples = samples,
                                   modifications = modifications)

  res_max_abn <- max(protein_map$munged_residues$abundance)

  if(!is.null(ncol) && !is.null(nrow)) {cli::cli_abort(c("x" = "only ncol or nrow can have a value, not both"))}

  p_map_all <- protein_map$munged_residues %>%
    ggplot2::ggplot() +
    ggplot2::geom_segment(ggplot2::aes(x=col,
                                       xend=col,
                                       y=row,
                                       yend=(row) - ((row_length * .67) * (abundance/res_max_abn))),
                          color=color_sequence, size=3.5)

  if(protein_map$munged_modifications %>% nrow() > 0) {
    p_map_all <- p_map_all +
      ggplot2::scale_color_manual(values = color_modifications) +
      ggplot2::geom_segment(data = protein_map$munged_modifications,
                            ggplot2::aes(x=col,
                                         xend=col,
                                         y=row,
                                         yend=(row) - ((row_length * .67) * (abundance/res_max_abn)),
                                         color=modification), size=3.5)

    if(show_modification_precent == TRUE){
      p_map_all <- p_map_all +
        ggplot2::geom_text(data = protein_map$munged_modifications,
                           ggplot2::aes(x=col,
                                        y=(row) - ((row_length * .67) * (abundance/res_max_abn)),
                                        label=paste0(signif(frequency*100,2),"%"),
                                        color=modification), hjust=0.33, vjust=-0.5, size=3)
    }
  }

  p_map_all <- p_map_all +
    ggplot2::geom_text(ggplot2::aes(col, row, label=aa), color='black', vjust = 1) +
    ggplot2::scale_x_continuous(breaks = c(1,1:(max(protein_map$munged_residues$col)/10) * 10)) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      panel.grid.minor.y=ggplot2::element_blank(),
      panel.grid.major.y=ggplot2::element_blank(),
      axis.title = ggplot2::element_blank()
    ) +
    ggplot2::facet_wrap(~sample, ncol=ncol, nrow=nrow) +
    ggplot2::scale_y_reverse(breaks = protein_map$munged_residues$row %>% unique()) +
    ggplot2::labs(title = protein)

  return(p_map_all)

}
