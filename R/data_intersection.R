#' Create a data subset
#'
#' @description
#' `intersection()` is a specalized function for sub-setting quantitative data
#' from a tidyproteomics data-object based data overlapping between sample groups.
#'
#' @param data tidyproteomics data object
#' @param .include when exporting the "intersection" this is the set of proteins
#' contained within the intersection of these samples
#' @param .exclude when exporting the "intersection" this is the set of proteins
#' found in these samples to exclude
#'
#' @return a tibble
#' @export
#'
#' @examples
#' library(dplyr, warn.conflicts = FALSE)
#' library(tidyproteomics)
#'
#' # creates a subset of just the proteins found in 'control'
#' hela_proteins %>%
#'    subset(imputed == 0) %>%
#'    intersection(.include = c('control'), .exclude = c('knockdown'))
#'
intersection <- function(
    data = NULL,
    .include = NULL,
    .exclude = NULL
){

  # visible bindings
  sample_id <- NULL

  check_data(data)
  identifier <- data$identifier

  if(identifier != 'protein'){ cli::cli_abort("This function currently only works for `proteins`") }

  sub_identifiers <- data %>% intersect_venn(include = .include, exclude = .exclude)
  cli::cli_alert_info("Intersection subset containing {.emph {length(sub_identifiers)}} {.emph {identifier}s}")

  for(x in c("quantitative", "accounting", "annotations")){
    w <- which(unlist(data[[x]][identifier]) %in% sub_identifiers)
    if(length(w) > 0){
      data[[x]] <- data[[x]][w,]
    }
  }

  w <- which(data$experiments$sample_id %in% unique(data$accounting$sample_id))
  if(length(w) > 0)
    data$experiments <- data$experiments[w,]

  return(data)
}
