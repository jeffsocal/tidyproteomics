#' Tidy-Quant data object print definition
#'
#' @param obj tidyproteomics data object
#'
#' @return print object summary
#'
#' @exportS3Method
#'
print.tidyproteomics <- function(
    obj
) {

  check_data(obj)
  obj_size <- as.numeric(object.size(obj))
  names_samples <- unique(unlist(obj$experiments$sample))
  names_accounting <- names(obj$accounting)
  names_accounting <- names_accounting[-which(grepl("sample|protein$|peptide$|modification$|group$", names_accounting))]

  dynamic_range <- obj %>% extract() %>%
    dplyr::group_by(sample) %>%
    dplyr::summarise(
      abundance_q005 = stats::quantile(abundance, .005, na.rm = TRUE),
      abundance_q995 = stats::quantile(abundance, .995, na.rm = TRUE),
      .groups = 'drop'
    )
  dynamic_range <- log10(mean(dynamic_range$abundance_q995, na.rm = TRUE)) - log10(mean(dynamic_range$abundance_q005, na.rm = TRUE))

  cli::cli_h2(cli::style_bold("{.emph Quantitative Proteomics Data Object}"))
  println("Origin", glue::glue("{obj$origin}"))
  println("", glue::glue("{obj$analyt} ({prettyunits::pretty_bytes(obj_size)})"))
  println("Quantitation", glue::glue("{nrow(obj$experiments)} files"))
  println("", glue::glue("{length(names_samples)} samples ({paste(names_samples, collapse=', ')})"))
  println("", glue::glue("{length(unique(unlist(obj$quantitative[obj$identifier[1]])))} {obj$identifier[1]}s"))
  println("", glue::glue("{signif(dynamic_range,2)} log10 dynamic range"))

  if(obj$quantitative_source != 'raw') {println(" *normalized", glue::glue("{obj$quantitative_source}"))}
  v_impt <- obj %>% get_variables('accounting')
  if('imputed' %in% v_impt) {
    v_impt <- obj$operations[which(grepl('imputed.*via', obj$operations))]
    v_impt <- sub(".*imputed\\s+", "", v_impt)
    v_impt <- gsub("(\\s|\\n)+", " ", v_impt)
    println(" *imputed", glue::glue("{v_impt}"))
  }

  println("Accounting", glue::glue("({length(names_accounting)}) {stringr::str_wrap(paste(names_accounting, collapse=' '), 60, exdent = 16)}"))

  if(!is.null(obj$annotations)) {
    names_annotations <- unique(unlist(obj$annotations$term))
    println("Annotations", glue::glue("({length(names_annotations)}) {stringr::str_wrap(paste(names_annotations, collapse=' '), 60, exdent = 16)}"))
  }

  if(!is.null(obj$analysis)) {
    println("Analyses", glue::glue("({length(obj$analysis)})"))
    for(comp in names(obj$analysis)) {
      analyses <- sort(names(obj$analysis[[comp]]), decreasing = TRUE)
      terms <- ''
      if('enrichment' %in% analyses) {terms <- glue::glue("({paste(names(obj$analysis[[comp]][['enrichment']]), collapse =', ')})")}
      println("", glue::glue("{comp} -> {paste(analyses, collapse=' & ')} {terms}"))
    }
  }

  println("")
  invisible(obj)

}


#' Tidy-Quant data object plot definition
#'
#' @param obj tidyproteomics data object
#'
#' @return print object summary
#'
#' @exportS3Method
#'
plot.tidyproteomics <- function(
    obj
) {
  plot_counts(obj)
  invisible(obj)
}

#' Helper function for printing messages
#'
#' @param name string
#' @param message string
#' @param pad_length string
#'
#' @return console print line
#'
println <- function(name = '',
                    message = '',
                    pad_length = 15) {
  cat(stringr::str_pad(name, pad_length, 'right'), message, "\n")
}
