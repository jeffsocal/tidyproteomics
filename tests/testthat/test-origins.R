test_that("origins assigns and updates sample_origin in experiments", {
  data(hela_proteins)

  # Column does not exist initially
  expect_false("sample_origin" %in% colnames(hela_proteins$experiments))

  # Assign single file origin
  res1 <- hela_proteins %>% origins(sample_file == "F1", .replace = "patient_1")
  expect_true("sample_origin" %in% colnames(res1$experiments))
  expect_equal(res1$experiments$sample_origin[res1$experiments$sample_file == "F1"], "patient_1")
  # Other runs should default to their sample_file
  expect_equal(res1$experiments$sample_origin[res1$experiments$sample_file == "F4"], "F4")

  # Assign multiple files by sample condition or regex
  res2 <- res1 %>% origins(sample == "control", .replace = "cohort_A")
  expect_true(all(res2$experiments$sample_origin[res2$experiments$sample == "control"] == "cohort_A"))
  expect_equal(res2$experiments$sample_origin[res2$experiments$sample_file == "F2"], "F2")

  # Operations recorded
  expect_true(any(grepl("sample_origin", unlist(res2$operations))))

  # Error on non-existent variable
  expect_error(hela_proteins %>% origins(non_existent == "foo", .replace = "bar"))

  # Error on non-matching value
  expect_error(hela_proteins %>% origins(sample_file == "F99", .replace = "bar"))
})
