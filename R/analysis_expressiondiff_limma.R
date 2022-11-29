#' Calculate expression differences between two-samples
#'
#' `expression_limma()` is a function for evaluating expression differences
#' between two sample sets via the limma algorithm
#'
#' @param data tidyproteomics data object
#' @param experiment a character string representing the experimental sample set
#' @param control a character string representing the control sample set
#'
#' @return a tibble
#'
expression_limma <- function(
    data = NULL,
    experiment = NULL,
    control = NULL
){

  # visible bindings
  identifier <- NULL
  samples <- NULL
  P.Value <- NULL
  adj.P.Val <- NULL
  logFC <- NULL
  AveExpr <- NULL
  B <- NULL
  log2_foldchange <- NULL
  average_expression <- NULL
  foldchange <- NULL
  proportional_expression <- NULL
  limma_t_statistic <- NULL
  abundance <- NULL
  abundance_log2 <- NULL
  imputed <- NULL
  n <- NULL

  data_quant <- data %>% extract(data$quantitative_source) %>% transform_log2()

  # only accept proteins with complete values
  l_comp_pro <- data_quant %>%
    dplyr::filter(abundance_log2 > 1) %>%
    dplyr::group_by(identifier) %>%
    dplyr::summarise(n = dplyr::n(),
                     .groups = 'drop') %>%
    dplyr::filter(n >= max(n) * .75) %>%
    dplyr::select(identifier) %>%
    unlist()

  data_quant_wide <- data_quant %>%
    dplyr::filter(sample %in% c(experiment, control)) %>%
    dplyr::filter(identifier %in% l_comp_pro) %>%
    tidyr::pivot_wider(names_from = c('sample', 'replicate'),
                       names_sep="___",
                       values_from = 'abundance_log2')

  data_quant_groups <- tibble::tibble(samples = colnames(data_quant_wide)[-1]) %>%
    tidyr::separate(samples, into = c('sample','replicate'), sep="\\_{3}", remove = F) %>%
    dplyr::group_by(sample) %>%
    dplyr::summarise(
      n_groups = length(sample), .groups = 'drop'
    )

  dm <- c()
  dm <- cbind(dm, rep(as.numeric(experiment == data_quant_groups$sample),
                      as.numeric(unlist(data_quant_groups[experiment = data_quant_groups$sample]$n_groups))))
  dm <- cbind(dm, rep(as.numeric(experiment != data_quant_groups$sample),
                      as.numeric(unlist(data_quant_groups[control = data_quant_groups$sample]$n_groups))))

  data_lmfit <- data_quant_wide %>%
    tibble::column_to_rownames('identifier') %>%
    as.data.frame() %>%
    limma::lmFit(dm) %>%
    limma::contrasts.fit(matrix(c(1,-1))) %>%
    limma::eBayes() %>%
    limma::topTable(number = 100000, coef = 1) %>%
    tibble::rownames_to_column("identifier") %>%
    tibble::as_tibble() %>%
    dplyr::rename(p_value = P.Value,
                  adj_p_value = adj.P.Val,
                  log2_foldchange = logFC,
                  average_expression = AveExpr,
                  limma_t_statistic = t,
                  limma_B_statistic = B) %>%
    dplyr::mutate(foldchange = invlog2(log2_foldchange),
                  average_expression = invlog2(average_expression)) %>%
    dplyr::mutate(proportional_expression = average_expression / sum(average_expression, na.rm = T)) %>%
    dplyr::relocate(foldchange, .after="log2_foldchange") %>%
    dplyr::relocate(proportional_expression, .after="average_expression") %>%
    dplyr::relocate(limma_t_statistic, .before="limma_B_statistic")

  return(data_lmfit %>% munge_identifier('separate', data$identifier))

}
