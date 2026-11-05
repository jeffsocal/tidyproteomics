test_that("subset handles various expressions including spaceless and negation", {
  data(hela_proteins)

  # Spaceless equality
  sub_ctl <- hela_proteins %>% subset(sample == "control", .verbose = FALSE)
  expect_true(all(sub_ctl$quantitative$sample == "control"))
  expect_equal(unique(sub_ctl$quantitative$sample), "control")

  # Standard string matching
  sub_ribo <- hela_proteins %>% subset(description %like% "Ribosome", .verbose = FALSE)
  expect_true(nrow(sub_ribo$quantitative) > 0)
  ribo_desc <- sub_ribo$annotations %>%
    dplyr::filter(term == "description") %>%
    dplyr::pull(annotation)
  expect_true(all(grepl("Ribosome", ribo_desc, ignore.case = TRUE)))

  # Negated string matching
  sub_noribo <- hela_proteins %>% subset(!description %like% "Ribosome", .verbose = FALSE)
  expect_true(nrow(sub_noribo$quantitative) > 0)
  noribo_desc <- sub_noribo$annotations %>%
    dplyr::filter(term == "description") %>%
    dplyr::pull(annotation)
  expect_false(any(grepl("Ribosome", noribo_desc, ignore.case = TRUE)))
})

test_that("subset on accounting keeps experiments and quantitative replicates synchronized", {
  data(hela_proteins)

  # Filter on accounting metric
  sub <- hela_proteins %>% subset(num_peptides >= 2, .verbose = FALSE)

  exp_map <- sub$experiments %>% dplyr::select(sample_id, sample, replicate) %>% dplyr::arrange(sample_id)
  quant_map <- sub$quantitative %>% dplyr::distinct(sample_id, sample, replicate) %>% dplyr::arrange(sample_id)

  expect_equal(exp_map$sample, quant_map$sample)
  expect_equal(exp_map$replicate, quant_map$replicate)

  # Verify downstream impute does not duplicate rows or produce key clashes
  n_rows_before <- nrow(sub$quantitative)
  imp <- sub %>% impute(base::min)
  expect_equal(nrow(imp$quantitative), n_rows_before)
})

