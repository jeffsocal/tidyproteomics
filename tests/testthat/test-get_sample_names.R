test_that("get_sample_names returns table by default and vector when requested", {
  data(hela_proteins)

  tbl_samples <- get_sample_names(hela_proteins)
  expect_s3_class(tbl_samples, "tbl_df")
  expect_true(all(c("sample", "replicates") %in% colnames(tbl_samples)))
  expect_equal(tbl_samples$sample, unique(hela_proteins$experiments$sample))

  vec_samples <- get_sample_names(hela_proteins, format = "vector")
  expect_true(is.character(vec_samples))
  expect_true(length(vec_samples) > 0)
  expect_true("control" %in% vec_samples || any(grepl("control", vec_samples)))
  expect_equal(vec_samples, unique(hela_proteins$experiments$sample))
})

