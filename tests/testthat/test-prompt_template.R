test_that("prompt_template returns valid prompt templates for all types", {
  # Test default (discovery)
  p_disc <- prompt_template("discovery")
  expect_true(is.character(p_disc))
  expect_true(grepl("What is biologically and technically interesting", p_disc))
  expect_true(grepl("tidyproteomics", p_disc))
  expect_true(grepl("SKILL\\.md", p_disc))
  expect_true(grepl("select_normalization", p_disc))
  expect_true(grepl("normalization check", p_disc, ignore.case = TRUE))
  expect_true(grepl("destination = \"png\"", p_disc))
  expect_true(grepl("without `ggsave\\(\\)`", p_disc))
  expect_true(grepl("plot_volcano", p_disc))

  # Test profiling
  p_prof <- prompt_template("profiling")
  expect_true(is.character(p_prof))
  expect_true(grepl("print\\(data\\)", p_prof))
  expect_true(grepl("data quality check", p_prof, ignore.case = TRUE))
  expect_true(grepl("summary\\(\"sample\"\\)", p_prof))
  expect_true(grepl("get_sample_names\\(\\)", p_prof))
  expect_true(grepl("show_annotations\\(\\)", p_prof))
  expect_true(grepl("summary\\(contamination = \"CRAP\"\\)", p_prof))

  # Test pairwise
  p_pair <- prompt_template("pairwise")
  expect_true(is.character(p_pair))
  expect_true(grepl("expression", p_pair))
  expect_true(grepl("TREATMENT_GROUP", p_pair))

  # Test multicontrast
  p_multi <- prompt_template("multicontrast")
  expect_true(is.character(p_multi))
  expect_true(grepl("analyze_expressions", p_multi))
  expect_true(grepl("analyze_enrichments", p_multi))

  # Test extracting
  p_ext <- prompt_template("extracting")
  expect_true(is.character(p_ext))
  expect_true(grepl("intersection", p_ext))
  expect_true(grepl("extract\\(na\\.rm = TRUE\\)", p_ext))
  expect_true(grepl("imputed == 0", p_ext))
  expect_true(grepl("plot_venn", p_ext))

  # Test invalid type argument
  expect_error(prompt_template("invalid_type"))
})
