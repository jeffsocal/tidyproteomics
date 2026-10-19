test_that("wilcoxon enrichment works properly and quickly", {
  d <- hela_proteins %>%
    expression(knockdown/control) %>%
    enrichment(knockdown/control, .terms = "biological_process", .method = "wilcoxon")

  res <- d$analysis[["knockdown/control"]]$enrichment[["biological_process"]]$data
  expect_s3_class(res, "tbl_df")
  expect_true(all(c("annotation", "p_value", "adj_p_value", "enrichment", "size") %in% colnames(res)))
  expect_true(nrow(res) > 0)
  expect_false("other" %in% res$annotation)
  expect_true(all(res$size >= 3))
  expect_true(all(res$p_value >= 0 & res$p_value <= 1))
  expect_true(all(res$adj_p_value >= 0 & res$adj_p_value <= 1))
})

test_that("wilcoxon handles empty annotation sets gracefully", {
  d <- hela_proteins %>% expression(knockdown/control)
  res <- tidyproteomics:::enrichment_wilcoxon(
    data_expression = d$analysis[["knockdown/control"]]$expression[0, ],
    data = d,
    term_group = "biological_process"
  )
  expect_s3_class(res, "tbl_df")
  expect_equal(nrow(res), 0)
  expect_true(all(c("annotation", "p_value", "adj_p_value", "enrichment", "size") %in% colnames(res)))
})
