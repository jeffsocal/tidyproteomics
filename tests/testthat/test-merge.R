test_that("merge works with base data frames and tidyproteomics objects", {
  # Base R merge compatibility
  df1 <- data.frame(id = 1:3, val1 = c("a", "b", "c"))
  df2 <- data.frame(id = 2:4, val2 = c("x", "y", "z"))
  merged_df <- merge(df1, df2, by = "id")
  expect_equal(nrow(merged_df), 2)
  expect_true(all(c("id", "val1", "val2") %in% colnames(merged_df)))

  # Tidyproteomics merge via S3 generic dispatch and function
  data(hela_proteins)
  p1 <- hela_proteins %>% subset(sample == "control")
  p2 <- hela_proteins %>% subset(sample == "knockdown")

  # S3 dispatch
  merged_tp <- merge(p1, p2)
  expect_true(inherits(merged_tp, "tidyproteomics"))
  expect_equal(merged_tp$origin, "Merged")
  expect_equal(length(unique(merged_tp$experiments$sample)), 2)

  # Explicit function
  merged_tp2 <- merge_tidyproteomics(list(p1, p2))
  expect_true(inherits(merged_tp2, "tidyproteomics"))
  expect_equal(merged_tp2$origin, "Merged")
})
