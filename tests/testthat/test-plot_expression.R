test_that("plot_expression creates expected ggplot and options work", {
  library(dplyr, warn.conflicts = FALSE)
  data(hela_proteins)

  hp_exp <- hela_proteins %>% expression(knockdown/control)

  # Basic plot with defaults (labels = protein)
  p1 <- plot_expression(
    hp_exp,
    knockdown/control,
    log2fc_min = 0.5,
    significance_column = "p_value",
    top_n = 4
  )
  expect_s3_class(p1, "ggplot")
  expect_equal(p1$labels$x, "Sample")
  expect_true(grepl("Knockdown / Control", p1$labels$title))

  # Faceting by gene_name
  p2 <- plot_expression(
    hp_exp,
    knockdown/control,
    log2fc_min = 0.5,
    significance_column = "p_value",
    top_n = 4,
    labels_column = "gene_name"
  )
  expect_s3_class(p2, "ggplot")

  # Regulation up / down / both
  p_up <- plot_expression(
    hp_exp,
    knockdown/control,
    log2fc_min = 0.5,
    significance_column = "p_value",
    top_n = 2,
    regulation = "up"
  )
  expect_s3_class(p_up, "ggplot")
  expect_true(grepl("\\(up\\)", p_up$labels$subtitle))

  p_down <- plot_expression(
    hp_exp,
    knockdown/control,
    log2fc_min = 0.5,
    significance_column = "p_value",
    top_n = 2,
    regulation = "down"
  )
  expect_s3_class(p_down, "ggplot")
  expect_true(grepl("\\(down\\)", p_down$labels$subtitle))

  # Fallback when no values meet threshold
  expect_message(
    p_fallback <- plot_expression(
      hp_exp,
      knockdown/control,
      log2fc_min = 20,
      significance_max = 1e-10,
      top_n = 2
    ),
    "no values significant"
  )
  expect_s3_class(p_fallback, "ggplot")
})
