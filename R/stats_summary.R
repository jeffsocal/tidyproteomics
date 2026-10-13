#' Summarize the protein accounting
#'
#' @description
#' `stats_summary()` is an analysis function that computes the protein summary
#' statistics for a given tidyproteomics data object.
#'
#' @param data tidyproteomics data object
#' @param group_by what to summarize
#'
#' @return a tibble
#'
stats_summary <- function(
    data,
    group_by = c("global", "sample", "replicate", "experiment", "normalization")
){

  # visible bindings
  sample_id <- NULL
  files <- NULL
  peptide <- NULL
  abundance <- NULL
  abundance_log2 <- NULL
  matches <- NULL
  identifier <- NULL
  sd <- NULL
  median <- NULL
  cv <- NULL
  CV <- NULL
  PC1 <- NULL
  PC2 <- NULL
  num_peptides <- NULL
  num_unique_peptides <- NULL
  peptides <- NULL
  quantifiable <- NULL
  protein <- NULL
  abundance_raw <- NULL
  protein_group <- NULL
  proteins <- NULL

  check_data(data)

  group_by <- rlang::arg_match(group_by)
  if(group_by == 'normalization') {
    quant_methods <- get_quant_names(data)
    tbl_stats <- lapply(quant_methods, function(qm) {
      d <- data
      d$quantitative_source <- qm
      res <- stats_summary(d, group_by = "global")
      dplyr::bind_cols(tibble::tibble(normalization = qm), res)
    }) %>% dplyr::bind_rows()
    return(tbl_stats)
  }

  if(group_by == 'experiment') { group_by <- 'sample_id' }
  if(group_by == 'global') { group_by <- NULL }
  identifier <- data$identifier

  cli::cli_progress_cleanup()
  cli::cli_progress_bar(type = 'tasks')
  if(data$analyte == 'peptides'){
    # cli::cli_progress_step(" ... calculating peptide -to- protein stats")
    tbl_long <- data %>% meld(single_quant_source = TRUE) %>%
      dplyr::rename(identifier = protein) %>%
      dplyr::group_by_at(c('sample_id', group_by, 'identifier')) %>%
      dplyr::summarise(
        num_peptides = dplyr::n(),
        num_unique_peptides = length(unique(peptide)),
        abundance = sum(abundance, na.rm = T),
        .groups = 'drop'
      )

    identifier <- 'protein'

  } else {
    tbl_long <- data %>% meld(single_quant_source = TRUE) %>%
      munge_identifier('combine', identifiers = identifier)

  }

  tbl_stats <- tbl_long %>%
    dplyr::group_by_at(group_by) %>%
    dplyr::summarise(
      identifiers = identifier %>% unique() %>% length(),
      peptides = ceiling(num_peptides %>% sum(na.rm = T) / length(unique(sample_id))),
      peptides_unique = ceiling(num_unique_peptides %>% sum(na.rm = T) / length(unique(sample_id))),
      quantifiable = signif(num_peptides[!is.na(abundance)] %>% sum(na.rm=T) / length(unique(sample_id)) / peptides, 3),
      .groups = 'drop'
    )

  tbl_cvs <- tbl_long %>%
    dplyr::group_by_at(c(group_by, 'identifier')) %>%
    dplyr::summarise(
      CV = signif(sd(abundance, na.rm = TRUE)/mean(abundance, na.rm = TRUE), 2),
      .groups = 'drop'
    ) %>%
    dplyr::group_by_at(group_by) %>%
    dplyr::summarise(
      CV = median(CV, na.rm = TRUE),
      .groups = 'drop'
    ) %>%
    dplyr::filter(!is.na(CV))

  calc_pca <- function(df) {
    tryCatch({
      df_log <- df %>%
        transform_log2(values = "abundance")

      mat <- df_log %>%
        dplyr::select(identifier, sample_id, abundance_log2) %>%
        tidyr::pivot_wider(
          names_from = "identifier",
          values_from = "abundance_log2",
          values_fill = 0
        ) %>%
        tibble::column_to_rownames("sample_id") %>%
        as.data.frame()

      mat <- mat[, which(apply(mat, 2, var, na.rm = TRUE) != 0), drop = FALSE]

      if (nrow(mat) < 2 || ncol(mat) < 2) {
        return(tibble::tibble(PC1 = NA_real_, PC2 = NA_real_))
      }

      pca <- stats::prcomp(mat, scale. = TRUE)
      imp <- base::summary(pca)$importance
      prop_var <- imp[which(rownames(imp) == "Proportion of Variance"), , drop = FALSE]

      pc1_val <- if ("PC1" %in% colnames(prop_var)) signif(prop_var[1, "PC1"], 2) else NA_real_
      pc2_val <- if ("PC2" %in% colnames(prop_var)) signif(prop_var[1, "PC2"], 2) else NA_real_

      tibble::tibble(PC1 = unname(pc1_val), PC2 = unname(pc2_val))
    }, error = function(e) {
      tibble::tibble(PC1 = NA_real_, PC2 = NA_real_)
    })
  }

  if (is.null(group_by)) {
    tbl_pca <- calc_pca(tbl_long)
    if (nrow(tbl_cvs) != 0) {
      tbl_stats <- dplyr::bind_cols(tbl_stats, tbl_cvs)
    }
    tbl_stats <- dplyr::bind_cols(tbl_stats, tbl_pca)
  } else {
    tbl_pca <- tbl_long %>%
      dplyr::group_by_at(group_by) %>%
      dplyr::group_modify(~ calc_pca(.x)) %>%
      dplyr::ungroup()

    if (nrow(tbl_cvs) != 0) {
      tbl_stats <- dplyr::left_join(tbl_stats, tbl_cvs, by = group_by)
    }
    tbl_stats <- dplyr::left_join(tbl_stats, tbl_pca, by = group_by)
  }

  w_cols <- which(grepl('^identifier[s_]*', colnames(tbl_stats)))
  colnames(tbl_stats)[w_cols] <- sub('identifier', identifier, colnames(tbl_stats)[w_cols])

  cli::cli_process_done()
  cli::cli_progress_cleanup()

  return(tbl_stats)
}
