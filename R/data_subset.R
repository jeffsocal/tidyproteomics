#' Create a data subset
#'
#' @description
#' `subset()` is the main function for sub-setting quantitative data from a tidyproteomics
#' data-object based on a regular expression and targeted annotation. This function
#' will return a smaller tidyproteomics data-object.
#'
#' Note: `rm.mbr()` is run as default, this is to remove MBR proteins that may no
#' longer have the original "anchor" observation present.
#'
#' @param data tidyproteomics data object
#' @param ... a three part expression (eg. x == a)
#' @param rm.mbr a boolean
#' @param .verbose a boolean
#'
#' @return a tibble
#' @export
#'
#' @examples
#' library(dplyr, warn.conflicts = FALSE)
#' library(tidyproteomics)
#'
#' # creates a subset of just Ribosomes, based on the string in the annotation
#' # protein_description
#' hela_proteins %>%
#'    subset(description %like% "Ribosome") %>%
#'    summary()
#'
#' # creates a subset without Ribosomes
#' hela_proteins %>%
#'    subset(!description %like% "Ribosome") %>%
#'    summary()
#'
subset <- function(
    data = NULL,
    ...,
    rm.mbr = TRUE,
    .verbose = TRUE
){

  # visible bindings
  sample_id <- NULL

  check_data(data)
  str_quo <- tidyproteomics_quo(...)
  if(is.null(str_quo)) { return(data) }

  variable <- str_quo[['variable']]
  operator <- str_quo[['operator']]
  value <- str_quo[['value']]
  inverse <- str_quo[['inverse']]
  inverse_str <- ''
  if(inverse == TRUE) { inverse_str <- '!' }

  identifier <- data$identifier
  which_segment <- get_segment(data, variable, .verbose)
  if(is.null(which_segment)) {return(data)}
  operation <- glue::glue("Data subset {inverse_str}`{variable}` {operator} `{value}`")

  if(.verbose == TRUE) {
    cli::cli_text("")
    cli::cli_div(theme = list(span.emph = list(color = "#ff4500"), span.info = list(color = "blue")))
    cli::cli_process_start("Subsetting data: {.emph {inverse_str}{variable}} {.info {operator}} {.emph {value}}")
  }

  identifier <- data$identifier

  if(which_segment == 'annotations') {

    select_by <- c(identifier,"sample_id")

    data[[which_segment]] <- data[[which_segment]] %>%
      tidyr::pivot_wider(id_cols = identifier,
                         names_from = 'term',
                         values_from = 'annotation') %>%
      down_select(str_quo) %>%
      tidyr::pivot_longer(cols = !dplyr::matches(paste0(identifier, '$')),
                          names_to = 'term', values_to = 'annotation')

    data$quantitative <- data$annotations %>%
      dplyr::select(identifier) %>% unique() %>%
      dplyr::inner_join(data$quantitative, by = identifier)

    data$experiments <- data$experiments %>%
      dplyr::filter(sample_id %in% data$quantitative$sample_id)

    data$accounting <- data$quantitative %>%
      dplyr::select(select_by) %>%
      dplyr::inner_join(data$accounting, by = select_by)

  } else {
    data[[which_segment]] <- down_select(data[[which_segment]], str_quo)
  }

  if(which_segment %in% c('experiments', 'quantitative')) {
    data$experiments <- data$experiments %>%
      dplyr::group_by(sample) %>%
      dplyr::mutate(replicate = as.character(dplyr::row_number())) %>%
      dplyr::ungroup()

    data$quantitative <- data$experiments %>%
      dplyr::select(c("sample_id","sample","replicate")) %>%
      dplyr::inner_join(data$quantitative %>% dplyr::select(!c('sample','replicate')),
                        by = "sample_id")

    data$accounting <- data$quantitative %>%
      dplyr::select(c("sample_id",identifier)) %>%
      dplyr::inner_join(data$accounting, by = c("sample_id",identifier))

  }

  if(which_segment == 'accounting') {

    data$experiments <- data$experiments %>%
      dplyr::filter(sample_id %in% unique(data$accounting$sample_id)) %>%
      dplyr::group_by(sample) %>%
      dplyr::mutate(replicate = as.character(dplyr::row_number())) %>%
      dplyr::ungroup()

    data$quantitative <- data$quantitative %>%
      dplyr::inner_join(data$accounting %>%
                          dplyr::select(c('sample_id',identifier)),
                        by = c('sample_id',identifier))
  }

  if(which_segment != 'accounting' & !is.null(data$annotations)) {
    data$annotations <- data$quantitative %>%
      dplyr::select(identifier) %>% unique() %>%
      dplyr::inner_join(data$annotations, by = identifier)
  }

  data$operations <- append(data$operations, operation)

  if(rm.mbr == TRUE) {data <- data %>% rm.mbr()}
  if(.verbose == TRUE) {cli::cli_process_done()}

  # cli::cli_alert_info("...[debug] {.emph {variable}} found in {.emph {which_segment}}")

  return(data)
}


#' Helper function to subset a data frame
#'
#' @param table a tibble
#' @param tidyproteomics_quo a character vector
#'
#' @return a tibble
#'
down_select <- function(
    table = NULL,
    tidyproteomics_quo = NULL
) {

  variable <- tidyproteomics_quo[['variable']]
  operator <- tidyproteomics_quo[['operator']]
  value <- tidyproteomics_quo[['value']]
  inverse <- tidyproteomics_quo[['inverse']]

  operator <- rlang::arg_match(operator, c("<","<=",">",">=","==","!=","%like%"))
  if(is.null(table) || is.null(variable) || is.null(value)) {return(table)}
  table <- as.data.frame(table)

  w <- c()
  if(operator == "<"){ w <- which(table[,variable] < value) }
  else if(operator == ">"){ w <- which(table[,variable] > value) }
  else if(operator == "<="){ w <- which(table[,variable] <= value) }
  else if(operator == ">="){ w <- which(table[,variable] >= value) }
  else if(operator == "=="){ w <- which(table[,variable] == value) }
  else if(operator == "!="){ w <- which(table[,variable] != value) }
  # if(operator == "%in%"){ w <- which(table[,variable] %in% value) }
  else if(operator == "%like%"){ w <- which(grepl(value, table[,variable], ignore.case = T)) }
  else {}

  if(length(w) > 0) {
    if(inverse == TRUE){ table <- table[-w,] }
    else{ table <- table[w,] }
  }

  return(tibble::as_tibble(table))
}


#' Helper function to subset a data frame
#'
#' @param ... a quo
#'
#' @return a list object
#'
tidyproteomics_quo <- function(...) {

  rlang_quo <- rlang::quo(...)
  quo_obj <- rlang::quo_text(rlang_quo)
  quo_obj <- sub("/", " / ", quo_obj)

  quo_obj <- sub("\\s+", " ", quo_obj)
  if(quo_obj == "") {return(NULL)}
  quo_str <- stringr::str_split(quo_obj, " ")[[1]]
  if(length(quo_str) < 3) {
    cli::cli_div(theme = list(span.emph = list(color = "#ff4500")))
    cli::cli_abort("improper expression from {.emph {quo_obj}}")
  }
  if(length(quo_str) > 3) {
    quo_str[3] <- paste(quo_str[3:length(quo_str)], collapse = " ")
    w <- 4:length(quo_str)
    quo_str <- quo_str[-w]
  }
  quo_str <- gsub('\"', '', quo_str)
  quo_str <- as.list(quo_str)
  quo_str[4] <- FALSE
  names(quo_str) <- c('variable','operator','value','inverse')

  # logic for "not"
  if(grepl("^\\!", quo_str[1])) {
    if(grepl("^\\%", quo_str[2])) { quo_str[4] <- TRUE }
    else { cli::cli_abort("undefined method for `!` at start of expression") }
    quo_str[1] <- sub("^\\!", "", quo_str[1])
  }

  # logic to turn characters back to a numeric
  if(grepl("[^0-9\\.]", quo_str[3])) {
    # keep as string
  } else { quo_str[3] <- as.numeric(quo_str[3]) }

  return(quo_str)
}

#' Helper function to get a name from the ...
#'
#' @param ... a quo
#'
#' @return a character string
#'
tidyproteomics_quo_name <- function(...){

  str_quo <- tidyproteomics_quo(...)
  if(is.null(str_quo)) { return(NULL) }
  return(paste(str_quo[['variable']], str_quo[['value']], sep="-"))
}

#' Helper function for subsetting
#'
#' @param a a dplyr tibble column reference
#' @param b a dplyr tibble column reference
#'
#' @return a character string
#'
#' @export
#'
`%like%` <- function(a, b) {
  grepl(b, a, ignore.case = T)
}

# #' @export
# `%!like%` <- function(a, b) {
#   !grepl(b, a, ignore.case = T)
# }
