test_that("stats_summary supports group_by = 'normalization' and reports CV, PC1, PC2", {
  data(hela_proteins)

  norm <- hela_proteins %>%
    normalize(.method = c("scaled", "median"))

  res <- tidyproteomics:::stats_summary(norm, group_by = "normalization")

  expect_s3_class(res, "tbl_df")
  expect_true("normalization" %in% colnames(res))
  expect_equal(res$normalization, get_quant_names(norm, format = "vector"))
  expect_true(all(c("proteins", "peptides", "peptides_unique", "quantifiable", "CV", "PC1", "PC2") %in% colnames(res)))
  expect_equal(nrow(res), length(get_quant_names(norm, format = "vector")))
  expect_true(all(!is.na(res$CV)))
  expect_true(all(!is.na(res$PC1)))
  expect_true(all(!is.na(res$PC2)))

  # Check column ordering: PC1 and PC2 should appear immediately after CV
  col_names <- colnames(res)
  cv_idx <- which(col_names == "CV")
  expect_equal(col_names[cv_idx + 1], "PC1")
  expect_equal(col_names[cv_idx + 2], "PC2")

  # Also test via summary(..., destination = 'return')
  res_summary <- norm %>% summary("normalization", destination = "return")
  expect_s3_class(res_summary, "tbl_df")
  expect_equal(res_summary$normalization, get_quant_names(norm, format = "vector"))
  expect_true(all(c("CV", "PC1", "PC2") %in% colnames(res_summary)))
})

test_that("stats_summary group_by = 'normalization' works on peptide data with PCA", {
  data(hela_peptides)

  norm_pep <- hela_peptides %>%
    normalize(.method = "scaled")

  res_pep <- tidyproteomics:::stats_summary(norm_pep, group_by = "normalization")
  expect_s3_class(res_pep, "tbl_df")
  expect_true("normalization" %in% colnames(res_pep))
  expect_equal(res_pep$normalization, get_quant_names(norm_pep, format = "vector"))
  expect_true(all(c("CV", "PC1", "PC2") %in% colnames(res_pep)))
  expect_true(all(!is.na(res_pep$PC1)))
})

test_that("stats_summary global and sample modes include CV, PC1, PC2", {
  data(hela_proteins)

  res_global <- tidyproteomics:::stats_summary(hela_proteins, group_by = "global")
  expect_true(all(c("CV", "PC1", "PC2") %in% colnames(res_global)))
  expect_false(is.na(res_global$PC1))

  res_sample <- tidyproteomics:::stats_summary(hela_proteins, group_by = "sample")
  expect_true(all(c("CV", "PC1", "PC2") %in% colnames(res_sample)))
  expect_equal(nrow(res_sample), 2)
  expect_false(any(is.na(res_sample$PC1)))
})
