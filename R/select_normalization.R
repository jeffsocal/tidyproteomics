#' Select a normalization method
#'
#' @description
#' `select_normalization()` selects the best normalization method base on low
#' CVs, low PCA (PC1), and wide Dynamic Range. This is a _passthrough_ function
#' as it returns the original tidyproteomics data-object.
#'
#' @param data tidyproteomics data object
#' @param normalization a character string
#'
#' @return a tidyproteomics data-object
#' @export
#'
#' @examples
#' library(dplyr, warn.conflicts = FALSE)
#' library(tidyproteomics)
#' hela_proteins <- hela_proteins %>%
#'   normalize(.method = c("scaled", "median", "linear", "limma", "loess","randomforest")) %>%
#'   select_normalization()
#'
select_normalization <- function(
    data = NULL,
    normalization = NULL
){

  # visible bindings
  origin <- NULL
  identifier <- NULL
  abundance_log2 <- NULL
  var <- NULL
  val <- NULL
  values <- NULL
  princom <- NULL
  abundance <- NULL
  abundance_95h <- NULL
  abundance_95l <- NULL
  abundance_sd <- NULL
  abundance_mean <- NULL
  abundance_cv <- NULL
  cv_mean <- NULL

  check_data(data)
  analyte <- data$analyte
  quant_values <- data %>% get_quant_names()

  if(!is.null(normalization)) {
    normalization <- rlang::arg_match(normalization, quant_values)
    data$quantitative_source <- normalization
    cli::cli_div(theme = list(span.emph = list(color = "blue")))
    cli::cli_alert_info("Normalization method set to {.emph {data$quantitative_source}}")
    data$operations <- append(data$operations, glue::glue("Normalization manually selected as {data$quantitative_source}."))
    return(data)
  }

  cli::cli_process_start("Selecting best normalization method")
  data_quant <- data %>% extract(values = quant_values)


  if(!c('origin') %in% names(data_quant)) {
    cli::cli_abort(c("x" = "Input {.emph data} has not been normalized to compare methods"))
  }

  tb_quant <- data_quant %>% transform_log2() %>%
    tidyr::unite(sample, sample, replicate)

  pca_sum <- list()
  for( value in quant_values ){
    df_pca <- tb_quant %>%
      dplyr::filter(origin == value) %>%
      dplyr::select(identifier, sample, abundance_log2) %>%
      tidyr::pivot_wider(
        names_from = "identifier",
        values_from = "abundance_log2",
        values_fill = 0
      ) %>%
      tibble::column_to_rownames('sample') %>%
      as.data.frame()

    df_pca <- df_pca[,which(apply(df_pca, 2, var, na.rm=TRUE) != 0)]

    pca_res <- stats::prcomp(df_pca, scale. = TRUE)
    pca_sum[[value]] <- base::summary(pca_res)$importance %>%
      tidyr::as_tibble() %>%
      dplyr::mutate(var = c('stdev','prop_var','cum_var')) %>%
      tidyr::pivot_longer(!var, names_to = 'princom', values_to = 'val') %>%
      dplyr::mutate(pcn = sub("[A-Z]+", "", princom) %>% as.numeric(),
                    values = value)
  }

  tb_final <- list(
    pca_sum %>%
      dplyr::bind_rows() %>%
      dplyr::filter(var == 'prop_var',
                    val != 0) %>%
      dplyr::group_by(values) %>%
      dplyr::mutate(val = cumsum(val)) %>%
      dplyr::ungroup() %>%
      dplyr::filter(princom == 'PC1') %>%
      dplyr::mutate(origin = values,
                    value = val,
                    metric = 'pca') %>%
      dplyr::mutate(value = min(value, na.rm = T) / value) %>%
      dplyr::select(c('origin','metric','value')),
    data_quant %>%
      dplyr::mutate(abundance = abundance %>% log10(),
                    origin = origin %>% as.factor()) %>%
      dplyr::filter(!is.infinite(abundance)) %>%
      dplyr::group_by(origin, sample) %>%
      dplyr::summarise(
        abundance_min = min(abundance, na.rm = T),
        abundance_max = max(abundance, na.rm = T),
        abundance_95l = stats::quantile(abundance, .025, na.rm=T),
        abundance_95h = stats::quantile(abundance, .975, na.rm=T),
        .groups = "drop"
      ) %>%
      dplyr::group_by(origin) %>%
      dplyr::summarise(
        value = mean(abundance_95h, na.rm = T) - mean(abundance_95l, na.rm = T),
        .groups = "drop"
      ) %>%
      dplyr::mutate(metric = 'dynamic_range',
                    value = value / max(value, na.rm = T)),
    data_quant %>%
      dplyr::mutate(origin = sub("data_vals_", "", origin),
                    origin = sub("\\.rds", "", origin),
                    origin = origin %>% as.factor()) %>%
      dplyr::group_by(origin, identifier, sample) %>%
      dplyr::summarise(
        abundance_mean = mean(abundance, na.rm = T),
        abundance_sd = stats::sd(abundance, na.rm = T),
        abundance_cv = abundance_sd / abundance_mean,
        .groups = "drop"
      ) %>%
      dplyr::ungroup() %>%
      dplyr::filter(!is.na(abundance_cv)) %>%
      dplyr::group_by(origin, sample) %>%
      dplyr::summarise(
        cv_mean = mean(abundance_cv, na.rm = T),
        cv_sd = stats::sd(abundance_cv, na.rm = T),
        .groups = "drop"
      ) %>%
      dplyr::group_by(origin) %>%
      dplyr::summarise(
        value = mean(cv_mean, na.rm = T),
        .groups = "drop"
      ) %>%
      dplyr::ungroup() %>%
      dplyr::filter(!is.na(value)) %>%
      dplyr::filter(!is.infinite(value)) %>%
      dplyr::mutate(metric = 'cvs',
                    value = min(value, na.rm = T) / value)) %>%
    dplyr::bind_rows() %>%
    dplyr::group_by(origin) %>%
    dplyr::summarise(value = sum(value, na.rm = T), .groups = 'drop') %>%
    dplyr::arrange(dplyr::desc(value))

  cli::cli_process_done()
  data$quantitative_source <- tb_final$origin[1]
  cli::cli_div(theme = list(span.emph = list(color = "blue")))
  cli::cli_alert_info(" ... selected {.emph {data$quantitative_source}}")

  data$operations <- append(data$operations, glue::glue("Normalization automatically selected as {data$quantitative_source}."))
  return(data)
}

