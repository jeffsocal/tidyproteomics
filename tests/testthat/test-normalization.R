test_that("normalization methods work correctly", {
  data(hela_proteins)

  # Test median normalization
  hp_med <- normalize(hela_proteins, .method = "median")
  expect_true(any(grepl("median", unlist(hp_med$operations))))
  expect_true("abundance_median" %in% colnames(hp_med$quantitative))

  # Test scaled normalization and verify linear sum invariance
  hp_scaled <- normalize(hela_proteins, .method = "scaled")
  expect_true("abundance_scaled" %in% colnames(hp_scaled$quantitative))
  sums <- hp_scaled$quantitative %>%
    dplyr::group_by(sample, replicate) %>%
    dplyr::summarise(tot = sum(abundance_scaled, na.rm = TRUE), .groups = "drop")
  expect_equal(length(unique(round(sums$tot, 2))), 1)
})
