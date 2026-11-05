test_that("plot_imputation_diagnostic generates expected ggplot and handles imputed features", {
  library(dplyr, warn.conflicts = FALSE)
  data(hela_proteins)

  hp <- hela_proteins %>%
    subset(description %like% "Ribosome|Proteasome|Histone") %>%
    impute(base::min, group_by_sample = TRUE) %>%
    expression(knockdown/control)

  # Test default automatic selection
  p1 <- plot_imputation_diagnostic(
    hp,
    knockdown/control,
    labels_column = "protein",
    destination = "plot",
    top_n = 4
  )
  expect_s3_class(p1, "ggplot")
  expect_true(grepl("Imputation Diagnostic", p1$labels$title))

  # Test nominated candidates
  cand <- head(unique(hp$quantitative[[hp$identifier[1]]]), 3)
  p2 <- plot_imputation_diagnostic(
    hp,
    knockdown/control,
    candidates = cand,
    labels_column = "protein",
    destination = "plot"
  )
  expect_s3_class(p2, "ggplot")
  expect_true(grepl("Imputation Diagnostic", p2$labels$title))
})
