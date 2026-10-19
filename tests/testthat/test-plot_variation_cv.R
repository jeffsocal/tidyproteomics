test_that("plot_variation_cv creates grob and saves cleanly without leaving Rplots.pdf", {
  data(hela_proteins)

  hp_norm <- hela_proteins %>% normalize(.method = "scaled")

  # Test destination = "plot"
  p <- plot_variation_cv(hp_norm, destination = "plot")
  expect_s3_class(p, "gtable")
  expect_s3_class(p, "grob")

  # Test destination = "png" does not create Rplots.pdf
  pdf_file <- "Rplots.pdf"
  png_file <- "plot_proteins_normalizated_variation.png"
  if (file.exists(pdf_file)) unlink(pdf_file)
  if (file.exists(png_file)) unlink(png_file)

  res <- plot_variation_cv(hp_norm, destination = "png")

  expect_s3_class(res, "tidyproteomics")
  expect_true(file.exists(png_file))
  expect_false(file.exists(pdf_file))

  # Clean up test artifact
  unlink(png_file)
  if (file.exists(pdf_file)) unlink(pdf_file)
})
