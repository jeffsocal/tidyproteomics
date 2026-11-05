test_that("plot_signatures generates expected ggplot and respects options", {
  library(dplyr, warn.conflicts = FALSE)
  library(tidyproteomics)
  data(hela_proteins)

  hp_exp <- hela_proteins %>% expression(knockdown/control)

  sigs <- list(
    "Ribosome" = c("RPL26", "RPS27A"),
    "Cell Cycle" = c("CKS2", "RACGAP1")
  )

  # Test basic lollipop plot with signatures
  p1 <- plot_signatures(
    hp_exp,
    knockdown/control,
    signatures = sigs,
    labels_column = "gene_name",
    plot_type = "lollipop",
    destination = "plot"
  )
  expect_s3_class(p1, "ggplot")
  expect_true(grepl("Biological Signatures", p1$labels$title))

  # Test bar plot
  p2 <- plot_signatures(
    hp_exp,
    knockdown/control,
    signatures = sigs,
    labels_column = "gene_name",
    plot_type = "bar",
    destination = "plot"
  )
  expect_s3_class(p2, "ggplot")

  # Test backward compatibility alias plot_programs
  p3 <- plot_programs(
    hp_exp,
    knockdown/control,
    programs = sigs,
    labels_column = "gene_name",
    destination = "plot"
  )
  expect_s3_class(p3, "ggplot")

  # Test set_a / set_b convenience interface
  p4 <- plot_signatures(
    hp_exp,
    knockdown/control,
    set_a = c("RPL26", "RPS27A"),
    set_b = c("CKS2", "RACGAP1"),
    label_a = "Ribosomal",
    label_b = "Growth",
    labels_column = "gene_name",
    destination = "plot"
  )
  expect_s3_class(p4, "ggplot")
})
