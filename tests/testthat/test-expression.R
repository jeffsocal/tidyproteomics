test_that("expression analysis and base compatibility work", {
  # Base R expression compatibility
  e <- expression(x + y)
  expect_true(is.expression(e))
  expect_equal(e, base::expression(x + y))

  # Tidyproteomics expression difference testing
  data(hela_proteins)
  hp_exp <- hela_proteins %>% expression(knockdown/control)
  expect_true("knockdown/control" %in% names(hp_exp$analysis))

  res <- hp_exp$analysis$`knockdown/control`$expression
  expect_true(all(c("log2_foldchange", "p_value", "adj_p_value") %in% colnames(res)))
  expect_true(all(res$p_value >= 0 & res$p_value <= 1, na.rm = TRUE))
  expect_true(all(res$adj_p_value >= 0 & res$adj_p_value <= 1, na.rm = TRUE))
})
