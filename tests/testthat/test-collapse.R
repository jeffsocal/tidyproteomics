test_that("collapse aggregates peptides to proteins with top_n options", {
  data(hela_peptides)

  # Default collapse (top_n = Inf)
  prot_all <- hela_peptides %>% collapse(.verbose = FALSE)
  expect_equal(prot_all$analyte, "proteins")
  expect_true(nrow(prot_all$quantitative) > 0)
  expect_true("protein" %in% colnames(prot_all$quantitative))

  # top_n = 2 collapse
  prot_top2 <- hela_peptides %>% collapse(top_n = 2, .verbose = FALSE)
  expect_equal(prot_top2$analyte, "proteins")
  expect_true(nrow(prot_top2$quantitative) > 0)
})
