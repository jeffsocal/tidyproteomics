#' A function for evaluating expression differences between two sample sets via the limma algorithm
#'
#' @param data tidyproteomics data object
#' @param experiment a character string representing the experimental sample set
#' @param control a character string representing the control sample set
#' @param .method a two-distribution test function returning a p_value for the null
#' hypothesis. Default is t.test. Example functions include t.test, wilcox.test,
#' stats::ks.test ...
#' @param .p.adjust a stats::p.adjust string for multiple test correction
#' @param ... pass through arguments
#'
#' @return a tibble
#'
expression_test <- function(
    data = NULL,
    experiment = NULL,
    control = NULL,
    .method = stats::t.test,
    ...,
    .p.adjust = 'BH'
){

  # visible bindings
  abundance <- NULL
  identifier <- NULL
  Var1 <- NULL
  Var2 <- NULL
  x <- NULL
  y <- NULL
  p_value <- NULL
  log2_foldchange <- NULL
  average_expression <- NULL
  proportional_expression <- NULL
  n_groups <- NULL
  foldchange <- NULL
  a_samples <- NULL
  imputed <- NULL

  check_data(data)
  experiment <- rlang::arg_match(experiment, unique(data$experiments$sample))
  control <- rlang::arg_match(control, unique(data$experiments$sample))

  if(control == experiment){
    cli::cli_div(theme = list(span.emph = list(color = "#ff4500"), span.info = list(color = "blue")))
    cli::cli_abort("x" = "Expression analysis must have different samples choosen",
                   "{.info experiment}:{.emph {experiment}} and {.info control}:{.emph {control}} are the same")
  }

  data_quant <- data %>% extract(data$quantitative_source) %>%
    dplyr::select(!dplyr::matches('^origin$'))

  data_account <- data$accounting %>%
    munge_identifier("combine", data$identifier) %>%
    dplyr::left_join(data$experiments, by = c('sample_id'))

  # only accept proteins with complete values
  l_comp_pro <- data_quant %>%
    dplyr::filter(abundance > 1) %>%
    dplyr::group_by(identifier) %>%
    dplyr::summarise(n_groups = dplyr::n(),
                     .groups = 'drop') %>%
    dplyr::filter(n_groups >= max(n_groups) * .75) %>%
    dplyr::select(identifier) %>%
    unlist()

  data_quant_wide <- data_quant %>%
    dplyr::filter(identifier %in% l_comp_pro) %>%
    dplyr::group_by(identifier, sample) %>%
    tidyr::nest() %>%
    tidyr::pivot_wider(names_from = 'sample', values_from = 'data')

  calc_fc <- function(x, y){
    y <- unlist(y$abundance)
    x <- unlist(x$abundance)
    y <- y[which(!is.na(y))]
    x <- x[which(!is.na(x))]

    if(length(x) == 0) return(NA)
    if(length(y) == 0) return(NA)

    tbl <- expand.grid(x, y) %>%
      dplyr::mutate(fc = Var1/Var2)

    return(log2(stats::median(tbl$fc)))
  }

  calc_pv <- function(x, y, .method = stats::t.test, ...){
    y <- unlist(y$abundance)
    x <- unlist(x$abundance)
    y <- y[which(!is.na(y))]
    x <- x[which(!is.na(x))]

    if(length(x) <= 1) return(1)
    if(length(y) <= 1) return(1)

    return(.method(y, x)$p.value)
  }

  calc_ae <- function(x, y){
    y <- c(unlist(y$abundance), unlist(x$abundance))
    y <- y[which(!is.na(y))]

    return(mean(y))
  }

  data_quant_out <- data_quant_wide %>%
    dplyr::rename(x = tidyselect::all_of(experiment)) %>%
    dplyr::rename(y = tidyselect::all_of(control)) %>%
    dplyr::mutate(log2_foldchange = purrr::map2(x, y, calc_fc),
                  p_value = purrr::map2(x, y, calc_pv, .method, ...),
                  average_expression = purrr::map2(x, y, calc_ae)) %>%
    tidyr::unnest(c(p_value, log2_foldchange, average_expression)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(foldchange = 2^(log2_foldchange),
                  proportional_expression = average_expression / sum(average_expression, na.rm=T)) %>%
    dplyr::select(identifier,
                  average_expression,
                  proportional_expression,
                  foldchange,
                  log2_foldchange,
                  p_value)

  data_quant_out$adj_p_value <- data_quant_out$p_value %>% stats::p.adjust(method = .p.adjust)

  return(data_quant_out %>% munge_identifier('separate', data$identifier))

}
