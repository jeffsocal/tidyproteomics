#' Retrieve standardized LLM analysis prompt templates
#'
#' @description
#' `prompt_template()` provides standardized, high-quality prompt templates designed
#' to be used alongside `SKILL.md` when querying Large Language Models (LLMs) such as
#' Claude, ChatGPT, Gemini, or local models (e.g. gpt-oss-20b, Ollama, LLaMA).
#' It facilitates a structured two-phase discovery workflow: Phase 1 (Reconnaissance & Profiling)
#' and Phase 2 (Biological Discovery: "What is interesting about this dataset?").
#'
#' @param type type of prompt template to retrieve:
#'   - `"discovery"`: Comprehensive "What is interesting about this dataset?" biological discovery prompt (Phase 2).
#'   - `"profiling"`: Dataset reconnaissance and metadata extraction prompt (Phase 1).
#'   - `"pairwise"`: Targeted two-sample differential abundance prompt.
#'   - `"multicontrast"`: Multi-condition contrast matrix / screen analysis prompt.
#'   - `"extracting"`: Venn diagram overlap and segment extraction prompt.
#'
#' @return invisible character string of the prompt template
#' @export
#'
#' @examples
#' library(tidyproteomics)
#' prompt_template("profiling")
#' prompt_template("discovery")
#' prompt_template("extracting")
#'
prompt_template <- function(
    type = c("discovery", "profiling", "pairwise", "multicontrast", "extracting")
){

  type <- rlang::arg_match(type)

  templates <- list(
    profiling = paste(
      'You are an expert computational proteomics specialist using the R package `tidyproteomics`.',
      'Refer to the instructions in SKILL.md.',
      '',
      'I have imported my quantitative mass spectrometry dataset in R, but need to inspect its structure before performing statistical analysis.',
      '',
      'Please provide a minimal, executable R script that:',
      '1. Loads `tidyproteomics` and `dplyr`.',
      '2. Imports the data file located at "[PATH_TO_DATA_FILE]" from platform "[PLATFORM: e.g. MaxQuant / ProteomeDiscoverer / FragPipe / DIA-NN]" as "[ANALYTE: proteins / peptides]".',
      '3. Runs `print(data)` to inspect the concise summary of data contents, sample names, number of proteins, imputed method, automatic normalization selection, and annotations.',
      '4. Runs `summary("sample")` as a data quality check to inspect replicate-level counts, CVs, and PCA separation.',
      '5. Prints all unique sample names using `get_sample_names()`.',
      '6. Prints available annotation categories using `show_annotations()`.',
      '7. Checks for potential sample contaminants using `summary(contamination = "CRAP")`.',
      '',
      'I will run this script in R and share the console output with you so we can design the comparative analysis.',
      sep = "\n"
    ),

    discovery = paste(
      'You are an expert quantitative proteomics data scientist. You have been provided with the `tidyproteomics` package instructions in SKILL.md.',
      '',
      'I have an imported `tidyproteomics` dataset with the following experimental setup:',
      '- Organism: [e.g. Homo sapiens / Mus musculus]',
      '- Biological Question: [e.g. Mechanism of action of Drug X / Knockout phenotype]',
      '- Sample Groups: [e.g. Treatment = "drug_10uM", Control = "dmso", 3 replicates each]',
      '- Available Annotations: [e.g. gene_name, description, biological_process, reactome_pathway]',
      '',
      'Your goal is to answer: "What is biologically and technically interesting about this dataset?"',
      '',
      'Generate a complete, self-contained R analysis script using `tidyproteomics` that:',
      '1. Evaluates Quality Control & Global Structure (Data Quality Check):',
      '   - Performs a data quality check with `summary("sample")` to inspect replicate CVs, quantifiable fraction, and PCA coordinates.',
      '   - Saves QC plots using built-in `destination = "png"` (e.g. `plot_variation_cv(destination = "png")`, `plot_dynamic_range(destination = "png")`, `plot_pca(destination = "png")`).',
      '   - Chains plot calls directly without `ggsave()`, capturing returned data objects in variables (`data <- data %>% plot_...`) to avoid console printing.',
      '2. Selects Optimal Normalization & Normalization Check:',
      '   - Tests candidate normalizations (`c("median", "linear", "loess", "limma")`) and selects the optimal method (`select_normalization()`).',
      '   - Runs `summary("normalization")` as a normalization check to confirm CV reduction and variance preservation across candidate normalizations.',
      '   - Verifies the automatic normalization selection via `print(data)`.',
      '3. Handles Missing Values:',
      '   - Imputes missing data appropriately (e.g. left-censored minimum within sample groups or random forest).',
      '   - Verifies the imputed method via `print(data)`.',
      '4. Uncovers Differential Expression:',
      '   - Compares Treatment vs Control using `expression(Treatment/Control, .method = "limma")`.',
      '   - Generates and saves a publication-ready volcano plot (`plot_volcano(destination = "png")`, log2fc_min = 1.0, significance_max = 0.05).',
      '5. Explores Biological Mechanism (Pathway Enrichment):',
      '   - Runs GSEA (`enrichment(Treatment/Control, .method = "gsea", .terms = "biological_process")`).',
      '   - Generates and saves an enrichment plot (`plot_enrichment(destination = "png")`).',
      '6. Highlights Key Findings:',
      '   - Extracts and prints the top 10 most significantly upregulated and downregulated proteins.',
      '   - Extracts and prints the top enriched biological pathways.',
      '   - Summarizes the findings into coherent biological hypotheses explaining what is most interesting about the dataset.',
      sep = "\n"
    ),

    pairwise = paste(
      'You are an expert quantitative proteomics data scientist using `tidyproteomics` (see SKILL.md).',
      '',
      'I need to analyze differential abundance between two specific sample groups:',
      '- Numerator (Experimental / Treatment): "[TREATMENT_GROUP]"',
      '- Denominator (Baseline / Control): "[CONTROL_GROUP]"',
      '- Statistical Method: [limma / t.test / wilcox.test]',
      '- P-value adjustment: [BH / bonferroni]',
      '',
      'Write a complete R script that:',
      '1. Normalizes the dataset and handles missing observations.',
      '2. Runs `expression([TREATMENT_GROUP]/[CONTROL_GROUP], .method = "[METHOD]")`.',
      '3. Plots a volcano plot highlighting top hits.',
      '4. Exports the differential results to an Excel spreadsheet (`export_analysis()`).',
      '5. Prints the audit trail with `operations()`.',
      sep = "\n"
    ),

    multicontrast = paste(
      'You are an expert quantitative proteomics data scientist using `tidyproteomics` (see SKILL.md).',
      '',
      'I have a multi-condition experiment with the following pairwise comparisons:',
      '[List comparisons, e.g.:',
      ' - c("knockdown", "control")',
      ' - c("rescue", "control")',
      ' - c("rescue", "knockdown")]',
      '',
      'Write a complete R script that:',
      '1. Normalizes the dataset using `select_normalization()`.',
      '2. Automates all two-sample tests in batch using `expression(.pairs = pairs, .method = "limma")`.',
      '3. Automates multi-term GSEA across biological processes and pathways using `enrichment(.pairs = pairs, .terms = c("biological_process", "reactome_pathway"))`.',
      '4. Saves out all figures and tables automatically using `analyze_expressions(labels_column = "gene_name")` and `analyze_enrichments(significance_max = 0.05)`.',
      '5. Generates a comparative expression plot (`plot_compexp()`) to contrast fold changes between conditions.',
      sep = "\n"
    ),

    extracting = paste(
      'You are an expert quantitative proteomics data scientist using `tidyproteomics` (see SKILL.md).',
      '',
      'I have an imported `tidyproteomics` dataset with multiple sample groups:',
      '[List groups, e.g. "control", "drug_low", "drug_high"]',
      '',
      'I want to investigate which proteins are uniquely expressed in specific conditions versus shared across groups.',
      '',
      'Write a complete R script that:',
      '1. Filters for non-imputed proteins using `subset(imputed == 0)` to ensure true mass spectrometry identifications.',
      '2. Generates and saves a Venn/Euler diagram (`plot_venn(destination = "png")` or `plot_euler(destination = "png")`).',
      '3. Uses `intersection(.include = ..., .exclude = ...)` to extract the underlying data for:',
      '   - Proteins unique to [GROUP_A] (absent in other groups).',
      '   - Proteins shared across all groups.',
      '   - Specific overlap segments of interest.',
      '4. Uses `extract(na.rm = TRUE)` to convert each extracted segment into a tidy long-format tibble and exports them to CSV (`readr::write_csv()`).',
      '5. Prints a concise summary of each extracted segment using `print()`.',
      sep = "\n"
    )
  )

  selected_prompt <- templates[[type]]

  cat(selected_prompt, "\n")
  invisible(selected_prompt)
}
