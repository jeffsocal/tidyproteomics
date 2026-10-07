#' Calculate expression differences between two-samples via MSstats mixed-effects model
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
expression_me <- function(
    data = NULL,
    experiment = NULL,
    control = NULL
){

  if (!requireNamespace("MSstats", quietly = TRUE)) {
    cli::cli_abort("Package {.pkg MSstats} is required for mixed-effects analysis. Please install it from Bioconductor.")
  }

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

  data_quant <- data %>% extract(data$quantitative_source) %>%
    dplyr::select(!dplyr::matches('^origin$'))

  # collect stats for downstream integration
  tbl_stats <- data %>% expression_test(experiment, control)

  # only accept proteins with complete values
  l_comp_pro <- data_quant %>%
    dplyr::filter(sample %in% c(experiment, control)) %>%
    dplyr::group_by(identifier, sample) %>%
    dplyr::summarise(n = dplyr::n(),
                     .groups = 'drop') %>%
    dplyr::group_by(identifier) %>%
    dplyr::summarise(min_group = min(n),
                     n = dplyr::n(),
                     .groups = 'drop') %>%
    dplyr::filter(n > 1, min_group > 0) %>%
    dplyr::select(identifier) %>%
    unlist()

  # inform of missing values if any
  if((length(unique(data_quant$identifier)) - length(l_comp_pro)) > 0){
    cli::cli_alert_warning("expression::mixed-effects removed {length(unique(data_quant$identifier)) - length(l_comp_pro)} proteins with completely missing values")
  }

  data_quant <- data_quant %>%
    dplyr::filter(sample %in% c(experiment, control)) %>%
    dplyr::filter(identifier %in% l_comp_pro) %>%
    dplyr::left_join(data$experiments, by = c("sample", "replicate"))

  if(!'sample_origin' %in% colnames(data_quant)){
    data_quant$sample_origin <- data_quant$sample_file
  }

  msstats_input <- data_quant %>%
    dplyr::transmute(
      ProteinName       = as.character(identifier),
      PeptideSequence   = as.character(identifier),      # Pseudo-peptide if starting from protein table
      PrecursorCharge   = 1,
      FragmentIon       = "NA",
      ProductCharge     = "NA",
      IsotopeLabelType  = "L",
      Condition         = as.character(sample),    # e.g., "Healthy", "Necrotic"
      BioReplicate      = as.character(sample_origin),   # CRITICAL: Maps multiple ROIs to Patient ID
      Run               = as.character(sample_file),       # Unique sample / LC-MS run name
      Fraction          = 1,
      # Convert log2 abundance to linear scale if tidyproteomics output is log2-transformed:
      Intensity         = ifelse(!is.na(abundance), abundance, NA_real_)
    )   %>%
    # Filter out missing/zero intensities if necessary
    dplyr::filter(!is.na(Intensity) & Intensity > 0)

  # ==============================================================================
  # 3. RUN MSSTATS DATA PROCESSING
  processed_data <- MSstats::dataProcess(
    raw                = msstats_input,
    normalization      = FALSE,   # Skip: already normalized via SVM in tidyproteomics
    summaryMethod      = "linear", # or "TMP" (Tukey's Median Polish)
    MBimpute           = FALSE,   # Skip: already imputed in tidyproteomics
    censoredInt        = "NA",
    logTrans           = 2,       # MSstats applies log2 transformation to linear intensities
    use_log_file       = FALSE
  )

  # ==============================================================================
  # 4. DEFINE CONTRAST MATRIX & RUN DIFFERENTIAL EXPRESSION
  unique_conditions <- levels(processed_data$ProteinLevelData$GROUP)

  contrast_matrix <- matrix(c(-1, 1), nrow = 1)
  colnames(contrast_matrix) <- c(as.character(experiment), as.character(control))
  rownames(contrast_matrix) <- as.character(glue::glue("{experiment}_vs_{control}"))

  # Run linear mixed-effects model comparison
  comparison_results <- MSstats::groupComparison(
    contrast.matrix = contrast_matrix,
    data            = processed_data,
    use_log_file    = FALSE
  )

  # ==============================================================================
  # 5. EXTRACT & INSPECT RESULTS TABLE
  diff_abundance_table <- comparison_results$ComparisonResult %>%
    dplyr::as_tibble() %>%
    dplyr::select(
      protein = Protein,
      adj_p_value = adj.pvalue,
      p_value = pvalue,
      log2_foldchange = log2FC,
      std_error = SE,
      test_stat = Tvalue,
      deg_freedom = DF,
      percent_missing = MissingPercentage
    ) %>%
    dplyr::filter(!is.infinite(log2_foldchange)) |>
    dplyr::left_join(tbl_stats |>
                       dplyr::select(matches(glue::glue("{paste(data$identifier, collapse='|')}|average_expression|proportional_expression"))),
                     by = data$identifier) |>
    dplyr::mutate(foldchange = invlog2(log2_foldchange)) %>%
    dplyr::arrange(adj_p_value) |>
    dplyr::mutate(proportional_expression = average_expression / sum(average_expression, na.rm = T)) %>%
    dplyr::relocate(foldchange, .after="log2_foldchange") %>%
    dplyr::relocate(proportional_expression, .after="average_expression") %>%

  return(diff_abundance_table %>% munge_identifier('separate', data$identifier))

}
