#' Create a data subset
#'
#' @description
#' `subset()` is the main function for sub-setting quantitative data from a tidyproteomics
#' data-object based on a regular expression and targeted annotation. This function
#' will return a smaller tidyproteomics data-object.
#'
#' @param data tidyproteomics data object
#' @param ... a three part expression (eg. x == a)
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
#'    subset(description %!like% "Ribosome") %>%
#'    summary()
#'
subset <- function(
    data = NULL,
    ...,
    .verbose = TRUE
){

  # visible bindings
  sample_id <- NULL

  str_quo <- tidyproteomics_quo(...)
  if(is.null(str_quo)) { return(data) }

  variable <- str_quo['variable']
  operator <- str_quo['operator']
  value <- str_quo['value']

  identifier <- data$identifier
  which_segment <- get_segment(data, variable, .verbose)
  if(is.null(which_segment)) {return(data)}
  operation <- glue::glue("Data subset `{variable}` {operator} `{value}`")

  if(.verbose == TRUE) {
    cli::cli_text("")
    cli::cli_div(theme = list(span.emph = list(color = "#ff4500"), span.info = list(color = "blue")))
    cli::cli_process_start("Subsetting data: {.emph {variable}} {.info {operator}} {.emph {value}}")
  }

  identifier <- data$identifier

  if(which_segment == 'annotations') {

    select_by <- c(identifier,"sample_id")

    data[[which_segment]] <- data[[which_segment]] %>%
      tidyr::pivot_wider(identifier, names_from = 'term', values_from = 'annotation') %>%
      down_select(variable, operator, value) %>%
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
    data[[which_segment]] <- down_select(data[[which_segment]], variable, operator, value)
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

  if(.verbose == TRUE) {cli::cli_process_done()}

  return(data)
}


#' Helper function to subset a data frame
#'
#' @param table a tibble
#' @param variable a character string
#' @param value a character string
#' @param operator a character string
#'
#' @return a tibble
#'
down_select <- function(
    table = NULL,
    variable = NULL,
    operator = c("<","<=",">",">=","==","!=","%like%","%!like%"),
    value = NULL
) {

  operator <- rlang::arg_match(operator)
  if(is.null(table) || is.null(variable) || is.null(value)) {return(table)}
  table <- as.data.frame(table)

  w <- c()
  if(operator == "<"){ w <- which(table[,variable] < value) }
  if(operator == ">"){ w <- which(table[,variable] > value) }
  if(operator == "<="){ w <- which(table[,variable] <= value) }
  if(operator == ">="){ w <- which(table[,variable] >= value) }
  if(operator == "=="){ w <- which(table[,variable] == value) }
  if(operator == "!="){ w <- which(table[,variable] != value) }
  if(operator == "%like%"){ w <- which(grepl(value, table[,variable], ignore.case = T)) }
  if(operator == "%!like%"){ w <- which(!grepl(value, table[,variable], ignore.case = T)) }

  if(length(w) > 0) {table <- table[w,]}

  return(tibble::as_tibble(table))
}


#' Helper function to subset a data frame
#'
#' @param table a tibble
#' @param variable a character string
#' @param value a character string
#' @param operator a character string
#'
#' @return a tibble
#'
tidyproteomics_quo <- function(...) {

  rlang_quo <- rlang::quo(...)
  quo_obj <- rlang::quo_text(rlang_quo)
  quo_obj <- sub("/", " / ", quo_obj)
  # if(rlang::is_empty(rlang::quo_get_env(rlang_quo)) == TRUE){
  #   # return(NULL)
  #   quo_obj <- sub("/", " / ", ...)
  # } else {
  #   quo_obj <- sub("/", " / ", rlang::quo_text(rlang_quo))
  # }
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
  names(quo_str) <- c('variable','operator','value')

  return(quo_str)
}

#' @export
`%like%` <- function(a, b) {
  grepl(b, a, ignore.case = T)
}

#' @export
`%!like%` <- function(a, b) {
  !grepl(b, a, ignore.case = T)
}
