test_that("get_sample_names returns unique sample names from data object", {
  data(hela_proteins)

  samples <- get_sample_names(hela_proteins)
  expect_true(is.character(samples))
  expect_true(length(samples) > 0)
  expect_true("control" %in% samples || any(grepl("control", samples)))
  expect_equal(samples, unique(hela_proteins$experiments$sample))
})
