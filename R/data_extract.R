#' Main function for extracting quantitative data from a tidyproteomics data-object
#'
#' @param data tidyproteomics data object
#' @param values character string vector
#' @param na.rm a boolean
#'
#' @return a tibble
#' @export
#'
extract <- function(
    data=NULL,
    values='raw',
    na.rm = FALSE
){

  # visible bindings
  abundance <- NULL
  origin <- NULL

  check_data(data)
  values <- rlang::arg_match(values, get_quant_names(data), multiple = T)
  if(!is.logical(na.rm) && !na.rm %in% c(TRUE, FALSE)){
    cli::cli_abort(c("x" = "na.rm must be TRUE or FALSE, not `{na.rm}`"))
  }

  d <- data$quantitative %>% dplyr::select(!"sample_id")
  d <- d %>% munge_identifier('combine', identifiers = data$identifier)

  if(values %>% length() == 1){
    d_out <- d %>%
      dplyr::select(!tidyselect::matches('abundance') | tidyselect::matches(paste0(values, "$"))) %>%
      dplyr::rename(abundance = tidyselect::matches(values)) %>%
      dplyr::mutate(origin = 'raw')
  } else {
    d_out <- d %>%
      tidyr::pivot_longer(
        cols = tidyselect::matches('abundance'),
        names_to = 'origin',
        values_to = 'abundance'
      ) %>%
      dplyr::mutate(origin = sub("abundance\\_", "", origin),
                    origin = origin %>% as.factor())

    d_out <- d_out %>%
      dplyr::mutate(origin = forcats::fct_relevel(origin, dplyr::intersect(values, names(d_out))))
  }

  if(na.rm == TRUE) {d_out <- d_out %>% dplyr::filter(!is.na(abundance), abundance > 0)}

  return(d_out)
}

#' Main function for munging peptide data from an extracted tidyproteomics data-object
#'
#' @param data tidyproteomics data object
#' @param munge character string vector (combine | separate)
#' @param identifiers a character vector of the identifiers
#'
#' @return a tibble
#'
munge_identifier <- function(
    data,
    munge = c("combine","separate"),
    identifiers = c('protein', 'peptide', 'modifications')
) {

  # visible bindings
  identifier <- NULL

  munge <- rlang::arg_match(munge)

  # unite or separate based on the identifier
  if( munge == 'combine'){
    data <- data %>% tidyr::unite(identifier,
                                  identifiers, sep=" | ")
  } else {
    data <- data %>% tidyr::separate(identifier,
                                     into = identifiers,
                                     sep="\\s\\|\\s",
                                     convert = TRUE)
  }
  return(data)
}
