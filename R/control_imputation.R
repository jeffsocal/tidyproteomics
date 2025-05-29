#' Main method for imputing missing values
#'
#' @param data a tidyproteomics list data-object
#' @param .function summary statistic function. Default is base::min, examples of
#' other functions include min, max, mean, sum. Note, NAs will be
#' be removed in the function call.
#' @param method a character string to indicate the imputation method (row, column, matrix).
#' Consider a data matrix of peptide/protein "rows" and dataset "columns". A 'row'
#' functions by imputing values between samples looking at the values for a given
#' peptide/protein, while the 'column' method imputes within a dataset of values.
#' The function 'randomforest' imputes using data from all rows and columns, or
#' the "matrix", without bias toward sample groups. If given a bias for
#' sample groups, expression differences would also bias sample groups. If it is the
#' case that sample groups should be biased (such as gene deletion), then it is
#' suggested to impute using min function and the 'within' method.
#' @param group_by_sample a boolean to indicate that the data should be grouped by
#' sample name to bias the imputation to within that sample.
#' @param cores the number of threads used to speed the calculation
#'
#' @return a tidyproteomics list data-object
#' @export
#'
#' @examples
#' library(dplyr, warn.conflicts = FALSE)
#' library(tidyproteomics)
#' hela_proteins %>% summary("sample")
#'
#' hela_proteins %>% impute(.function = stats::median) %>% summary("sample")
#'
#' hela_proteins %>% impute(.function = impute.randomforest, method = "matrix") %>% summary("sample")
#'
impute <- function(
    data = NULL,
    .function = base::min,
    method = c('row', 'column', 'matrix'),
    group_by_sample = FALSE,
    cores = 2
){

  # visible bindings
  col_num <- NULL
  abundance <- NULL
  imputed <- NULL
  imputed.y <- NULL

  method <- rlang::arg_match(method)
  if(!is.logical(group_by_sample)) { cli::cli_abort("group_by_sample must be one of [TRUE, FALSE]") }
  group_by_sample <- group_by_sample[1]

  check_data(data)
  quant_source <- data$quantitative_source
  identifier <- data$identifier

  cli::cli_progress_bar(type = 'tasks')

  make_matrix <- function(x){
    x %>% tidyr::pivot_wider(names_from = dplyr::matches('sample|replicate'),
                             values_from = 'abundance',
                             names_sep = "-sep-") %>%
      tibble::column_to_rownames("identifier") %>% as.matrix()
  }

  undo_matrix <- function(x, group_by_sample = FALSE){
    sample_replicate <- NULL
    if(group_by_sample == TRUE){
      x %>% as.data.frame() %>%
        tibble::rownames_to_column("identifier") %>%
        tidyr::pivot_longer(cols = !identifier,
                            names_to = 'replicate',
                            values_to = 'abundance')
    } else {
      x %>% as.data.frame() %>%
        tibble::rownames_to_column("identifier") %>%
        tidyr::pivot_longer(cols = !identifier,
                            names_to = 'sample_replicate',
                            values_to = 'abundance') %>%
        tidyr::separate(sample_replicate,
                        into = c('sample','replicate'),
                        sep = "-sep-")
    }
  }

  impute_function_str <- enquote(.function) %>% paste(collapse = ' ') %>% stringr::str_replace_all("\\s+", " ")
  cli::cli_div(theme = list(span.emph = list(color = "#ff4500")))
  cli::cli_progress_step("Imputing by {.emph {method}} using the function {.emph {impute_function_str}}")

  all_na <- function(x){ length(which(!is.na(x))) == 0 }

  # create a matrix for imputing,
  table <- data %>%
    extract(quant_source) %>%
    dplyr::select(!dplyr::matches("^origin$"))

  if(group_by_sample == TRUE){
    # group by sample ?
    table <- table %>%
      dplyr::group_by(sample) %>%
      tidyr::nest() %>%
      dplyr::mutate(data = purrr::map(data, make_matrix))
  } else {
    table <- table %>%
      tidyr::nest(data = tidyr::everything()) %>%
      dplyr::mutate(data = purrr::map(data, make_matrix))
  }

  for(i in 1:length(table$data)){
    this_matrix <- table$data[[i]]

    if(method == 'row'){
      this_matrix_na <- which(this_matrix %>% base::apply(1, all_na) == TRUE)
      if(length(this_matrix_na) > 0){ this_matrix <- this_matrix[-this_matrix_na,] }
    }

    # find the values to impute
    w_impute <- which(is.na(this_matrix), arr.ind = TRUE)
    if(length(w_impute) == 0) { next() }

    # impute by row or col
    # =================================
    if(method == 'row'){
      row_val <- this_matrix %>% base::apply(1, function(x, .f = .function){ .f(as.vector(stats::na.omit(x))) })
      this_matrix[w_impute] <- row_val[w_impute[,"row"]]
    }

    if(method == 'column'){
      col_val <- this_matrix %>% base::apply(2, function(x, .f = .function){ .f(as.vector(stats::na.omit(x))) })
      this_matrix[w_impute] <- col_val[w_impute[,"col"]]
    }
    # =================================

    if(method == 'matrix'){

      n_row <- nrow(this_matrix)
      n_col <- ncol(this_matrix)
      this_matrix <- .function(as.matrix(this_matrix), cores = cores)

      if(is.null(nrow(this_matrix)) | is.null(ncol(this_matrix))) {
        cli::cli_abort("Imputation resulted in a non-matrix - check your function to ensure it returns a matrix")
      }
      if(n_row != nrow(this_matrix) | n_col != ncol(this_matrix)){
        cli::cli_abort("Imputation resulted in an irregular matrix - check your function to ensure it returns a matrix of the same size")
      }

    }

    # account for the imputed values
    this_impute <- data.frame(
      identifier = rownames(w_impute),
      col_num = w_impute[,2],
      imputed = TRUE
    ) %>%
      dplyr::left_join(
        tibble::tibble(replicate = colnames(this_matrix)) %>%
          dplyr::mutate(col_num = 1:dplyr::n()),
        by = 'col_num'
      ) %>%
      dplyr::select(!col_num) %>%
      tidyr::pivot_wider(names_from = 'replicate', values_from = 'imputed') %>%
      as.data.frame() %>%
      tibble::column_to_rownames('identifier') %>%
      undo_matrix(group_by_sample) %>%
      dplyr::rename(imputed = abundance) %>%
      dplyr::filter(!is.na(imputed))

    this_matrix <- this_matrix %>%
      undo_matrix(group_by_sample)

    # join the values and the markers for imputation
    this_matrix <- this_matrix %>%
      dplyr::left_join(
        this_impute,
        by = setdiff(colnames(this_matrix), c('abundance')),
      ) %>%
      dplyr::mutate(imputed = ifelse(is.na(imputed), FALSE, imputed))

    table$data[[i]] <- this_matrix
  }

  table <- table %>%
    tidyr::unnest(data) %>%
    dplyr::ungroup()

  table <- data$experiments %>%
    dplyr::select(c('sample_id', 'sample', 'replicate')) %>%
    dplyr::full_join(table, by = c('sample', 'replicate'))

  data <- data %>%
    merge_quantitative(table %>% dplyr::select(!imputed), quant_source)

  table <- table %>%
    munge_identifier("separate", data$identifier) %>%
    dplyr::select(!dplyr::matches('abundance')) %>%
    dplyr::select(!dplyr::matches('sample_id')) %>%
    dplyr::full_join(data$experiments, by = c('sample', 'replicate')) %>%
    dplyr::select(!c('sample', 'replicate', 'import_file', 'sample_file'))

  join_names <- colnames(table)[-which("imputed" == colnames(table))]
  data$accounting <- data$accounting %>%
    dplyr::full_join(table, by = join_names)

  if('imputed.y' %in% colnames(data$accounting)){
    data$accounting <- data$accounting %>%
      dplyr::rename(imputed = imputed.y) %>%
      dplyr::select(!c('imputed.x'))
  }

  cli::cli_progress_bar(type = 'tasks')

  impute_msg <- glue::glue("... {data$accounting %>% dplyr::filter(imputed == TRUE) %>% nrow()} values imputed")
  cli::cli_alert_info(impute_msg)
  data$operations <- append(data$operations, glue::glue("Missing values imputed by '{method}' samples via '{impute_function_str}' group_by_sample '{group_by_sample}'."))
  data$operations <- append(data$operations, impute_msg)

  return(data)
}
