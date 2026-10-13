test_that("export_llm_skill exports SKILL.md and contains required sections", {
  temp_dir <- tempfile("test_skill_dir")
  dir.create(temp_dir)
  on.exit(unlink(temp_dir, recursive = TRUE), add = TRUE)

  # Test default export into temp_dir
  saved_path <- export_llm_skill(destination = temp_dir)
  expect_true(file.exists(saved_path))
  expect_equal(basename(saved_path), "SKILL.md")

  # Read exported content
  content <- readLines(saved_path, warn = FALSE)
  content_text <- paste(content, collapse = "\n")

  # Verify frontmatter and required sections
  expect_true(grepl("name:\\s*tidyproteomics", content_text))
  expect_true(grepl("supported_formats", content_text, ignore.case = TRUE))
  expect_true(grepl("api_reference", content_text, ignore.case = TRUE))
  expect_true(grepl("ProteomeDiscoverer", content_text))
  expect_true(grepl("MaxQuant", content_text))
  expect_true(grepl("export_llm_skill", content_text))

  # Test custom file name
  custom_path <- export_llm_skill(destination = temp_dir, file_name = "CUSTOM_SKILL.md")
  expect_true(file.exists(custom_path))
  expect_equal(basename(custom_path), "CUSTOM_SKILL.md")

  # Test overwrite = FALSE when file already exists
  expect_error(
    export_llm_skill(destination = temp_dir, file_name = "CUSTOM_SKILL.md", overwrite = FALSE),
    "already exists"
  )
})

test_that("path_to_package_data can resolve skill", {
  skill_path <- path_to_package_data("skill")
  expect_true(file.exists(skill_path))
  expect_true(grepl("SKILL\\.md$", skill_path, ignore.case = TRUE))
})
