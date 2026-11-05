test_that("plot_candidates generates expected ggplot and respects status tags", {
  library(dplyr, warn.conflicts = FALSE)
  library(tidyproteomics)
  data(hela_proteins)

  hp_exp <- hela_proteins %>% expression(knockdown/control)

  # Test nominated candidates with status mapping
  status_map <- c("RPL26" = "Validated", "RPS27A" = "Validated", "CKS2" = "Discordant")

  p1 <- plot_candidates(
    hp_exp,
    knockdown/control,
    candidates = names(status_map),
    status = status_map,
    labels_column = "gene_name",
    destination = "plot"
  )
  expect_s3_class(p1, "ggplot")
  expect_true(grepl("Candidate Targets", p1$labels$title))

  # Test backward compatibility alias plot_forest
  p2 <- plot_forest(
    hp_exp,
    knockdown/control,
    proteins = names(status_map),
    status = status_map,
    labels_column = "gene_name",
    destination = "plot"
  )
  expect_s3_class(p2, "ggplot")

  # Test default fallback top candidates
  p3 <- plot_candidates(
    hp_exp,
    knockdown/control,
    labels_column = "gene_name",
    destination = "plot"
  )
  expect_s3_class(p3, "ggplot")
})
