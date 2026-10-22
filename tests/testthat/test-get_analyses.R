test_that("get_analyses returns empty table/vector when no analyses present", {
  data(hela_proteins)

  tbl_empty <- get_analyses(hela_proteins)
  expect_s3_class(tbl_empty, "tbl_df")
  expect_equal(nrow(tbl_empty), 0)
  expect_equal(colnames(tbl_empty), c("comparison", "analysis", "term", "method"))

  vec_empty <- get_analyses(hela_proteins, format = "vector")
  expect_equal(vec_empty, character())
})

test_that("get_analyses returns correct table and vector when analyses exist", {
  data(hela_proteins)

  suppressMessages({
    d <- hela_proteins %>%
      expression(knockdown/control, .method = stats::t.test)
  })

  tbl_exp <- get_analyses(d)
  expect_s3_class(tbl_exp, "tbl_df")
  expect_equal(nrow(tbl_exp), 1)
  expect_equal(tbl_exp$comparison, "knockdown/control")
  expect_equal(tbl_exp$analysis, "expression")
  expect_true(is.na(tbl_exp$term))

  vec_exp <- get_analyses(d, format = "vector")
  expect_equal(vec_exp, "knockdown/control:expression")

  suppressMessages({
    d <- d %>%
      enrichment(knockdown/control, .terms = "biological_process", .method = "wilcoxon")
  })

  tbl_both <- get_analyses(d)
  expect_s3_class(tbl_both, "tbl_df")
  expect_equal(nrow(tbl_both), 2)
  expect_equal(tbl_both$comparison, c("knockdown/control", "knockdown/control"))
  expect_equal(tbl_both$analysis, c("expression", "enrichment"))
  expect_equal(tbl_both$term, c(NA_character_, "biological_process"))

  vec_both <- get_analyses(d, format = "vector")
  expect_equal(vec_both, c("knockdown/control:expression", "knockdown/control:enrichment:biological_process"))
})

test_that("all other get_* functions support table and vector formats", {
  data(hela_proteins)

  # get_quant_names
  tbl_quant <- get_quant_names(hela_proteins)
  expect_s3_class(tbl_quant, "tbl_df")
  expect_true(all(c("quant_name", "selected") %in% colnames(tbl_quant)))
  vec_quant <- get_quant_names(hela_proteins, format = "vector")
  expect_true(is.character(vec_quant))
  expect_equal(tbl_quant$quant_name, vec_quant)

  # get_accountings
  tbl_acc <- get_accountings(hela_proteins)
  expect_s3_class(tbl_acc, "tbl_df")
  expect_true(all(c("accounting", "data_type") %in% colnames(tbl_acc)))
  vec_acc <- get_accountings(hela_proteins, format = "vector")
  expect_true(is.character(vec_acc))
  expect_equal(tbl_acc$accounting, vec_acc)

  # get_annotation_terms
  tbl_terms <- get_annotation_terms(hela_proteins)
  expect_s3_class(tbl_terms, "tbl_df")
  expect_true(all(c("term", "type", "unique_annotations") %in% colnames(tbl_terms)))
  vec_terms <- get_annotation_terms(hela_proteins, format = "vector")
  expect_true(is.character(vec_terms))
  expect_equal(tbl_terms$term, vec_terms)

  # get_annotations
  tbl_anno <- get_annotations(hela_proteins, "biological_process")
  expect_s3_class(tbl_anno, "tbl_df")
  vec_anno <- get_annotations(hela_proteins, "biological_process", format = "vector")
  expect_true(is.character(vec_anno))
  expect_equal(sort(unique(tbl_anno$annotation)), sort(vec_anno))

  # get_variables
  tbl_vars <- get_variables(hela_proteins, "experiments")
  expect_s3_class(tbl_vars, "tbl_df")
  expect_true(all(c("segment", "variable") %in% colnames(tbl_vars)))
  vec_vars <- get_variables(hela_proteins, "experiments", format = "vector")
  expect_true(is.character(vec_vars))
  expect_equal(tbl_vars$variable, vec_vars)

  # get_unique_variables
  tbl_uvar <- get_unique_variables(hela_proteins, "sample")
  expect_s3_class(tbl_uvar, "tbl_df")
  expect_true("value" %in% colnames(tbl_uvar))
  vec_uvar <- get_unique_variables(hela_proteins, "sample", format = "vector")
  expect_true(is.character(vec_uvar))
  expect_equal(tbl_uvar$value, vec_uvar)
})
