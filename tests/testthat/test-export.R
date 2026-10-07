test_that("export_quant outputs valid tables and handles scaling", {
  data(hela_proteins)

  # Default export: cleanly exports sample names without redundant abundance_scaled-none
  eq_def <- hela_proteins %>% export_quant()
  expect_false(any(grepl("abundance_scaled-none", colnames(eq_def))))
  expect_true(any(grepl("control_1", colnames(eq_def))))

  # Scaled between export: produces non-NA values and scaled columns
  eq_btw <- hela_proteins %>% export_quant(scaled = "between")
  cols_scaled <- colnames(eq_btw)[grepl("abundance_scaled", colnames(eq_btw))]
  expect_true(length(cols_scaled) > 0)
  expect_false(all(is.na(eq_btw[[cols_scaled[1]]])))
})
