---
name: tidyproteomics
description: Comprehensive expert guide, complete API reference, supported formats, and canonical workflows for quantitative proteomics post-analysis in R using tidyproteomics. Use when importing mass spectrometry outputs (ProteomeDiscoverer, MaxQuant, FragPipe, DIA-NN, Skyline, mzTab), performing exploratory data analysis, filtering, normalization, imputation, differential expression analysis (limma, t-test, wilcoxon), pathway and GO term enrichment analysis (GSEA, Fisher's exact, Wilcoxon), batch automated pipelines (.pairs, .terms, analyze_expressions, analyze_enrichments), Venn/Euler segment extraction (intersection, extract), peptide-to-protein collapsing, sequence visualization, or formatting tidyproteomics S3 objects.
---

# tidyproteomics: LLM Agent & Developer Skill Guide

`tidyproteomics` is an R package that introduces an S3 data object architecture and tidyverse-aligned pipeable framework for the post-analysis of quantitative proteomics mass spectrometry data. It simplifies importing, filtering, normalizing, imputing, analyzing, and visualizing proteomics datasets while maintaining a complete, non-destructive audit trail of all operations.

---

## 1. Core Architecture & Philosophy

### The `tidyproteomics` S3 Object Structure
A `tidyproteomics` object is an S3 list containing:
- `$analyte`: Character string declaring the measurement level: `"proteins"` or `"peptides"`.
- `$quantitative`: Tibble storing quantitative values per observation across samples and replicates (e.g., `sample`, `replicate`, `identifier`, `abundance_raw`, `abundance_normalized`, etc.).
- `$annotations`: Tibble mapping identifiers to biological annotations (e.g., `protein`, `gene_name`, `description`, `biological_process`, `cellular_component`, `molecular_function`, `wiki_pathway`, `reactome_pathway`).
- `$accounting`: Tibble storing peptide/PSM counts, identification confidence, and match-between-runs (MBR) flags.
- `$experiments`: Tibble detailing sample metadata, condition labels, and raw data file associations.
- `$operations`: Character list logging every transformation applied to the dataset (audit trail).
- `$analysis`: List storing completed statistical comparisons, differential expression results, and enrichment computations.

### Key Principles
1. **Pipe-First Design**: All primary operations accept the `data` object as their first argument and return an updated `tidyproteomics` object, designed for magrittr `%>%` pipelines.
2. **Concise Object State via `print(data)`**: Calling `print(data)` (or entering `data` at the console) exposes a concise, structured summary displaying data contents (origin, analyte, file count, memory footprint, dynamic range, % missing values), sample names, number of proteins/peptides, imputed method (e.g., `*imputed via randomforest`), automatic normalization selection (e.g., `*normalized loess`), and biological annotations.
3. **Data Quality & Normalization Checks via `summary(data)`**: `summary(data)` is the core quality assessment tool. Use `summary(data)` or `summary(data, "sample")` as a data quality check (replicate counts, quantifiable fraction, CV, and PC1 & PC2 separation), `summary(data, contamination = 'CRAP')` as a contaminant quality check, and `summary(data, "normalization")` as a normalization check to compare candidate normalizations on CV reduction and PCA variance.
4. **Non-Destructive Transformations**: Normalizations and imputations do not overwrite raw values by default; they add new columns or states (e.g., `abundance_normalized`), tracked via `data$quantitative_source`.
5. **Audit Trail**: Every filtering, normalization, imputation, and analysis step is automatically appended to `$operations`. Retrieve this history at any time using `operations(data)`.
6. **Quasiquotation & Tidy Evaluation**: Functions like `subset()`, `expression()`, and `export_analysis()` support unquoted expressions, formulas, and tidyverse syntax.
7. **Plotting Pipeline & Built-In Saving (`plot_save`)**: All `plot_*()` functions inherit saving arguments directly from `plot_save()` (`destination = c("plot", "save", "png", "svg", "tiff", "jpeg")`, `height`, `width`). When saving (`destination = "save"`, `"png"`, etc.), the plot is saved to disk and the function returns the **`tidyproteomics` data object**, enabling pipeline concatenation of multiple plots at once (`data %>% plot_counts(destination = "png") %>% plot_pca(destination = "png")`). However, if not captured (assigned or wrapped in `invisible()`), the final returned `data` object will trigger `print(data)` at the console. **Never apply `ggsave()` to `plot_*()` outputs**—saving is handled internally, and saving mode returns the `tidyproteomics` data object, not a ggplot object.
8. **Batch Automation Architecture (`automating.Rmd`)**: Native support for batch differential analysis across arbitrary contrast matrices via `.pairs = list(c('A', 'control'), c('B', 'control'), ...)` and multi-term pathway annotations via `.terms = c('biological_process', 'reactome_pathway', ...)`. End-to-end pipelines seamlessly chain into `analyze_expressions()` and `analyze_enrichments()` to automatically output publication-ready volcano plots, multi-term summary bar charts, and CSV tables for all pairs directly to disk in a single execution.
9. **Visual-to-Analytical Data Extraction (`extracting.Rmd`)**: Transition from high-level Venn/Euler visualizations (`plot_venn()`, `plot_euler()`) to targeted data extraction of specific overlapping or exclusive segments via `intersection(.include = ..., .exclude = ...)`. Because `intersection()` returns a valid, filtered `tidyproteomics` object (with synchronized quantitative, accounting, and experiment files), you can immediately pipe it into downstream analyses or extract tidy long tibbles via `extract(na.rm = TRUE)`.


---

## 2. Supported Formats & Data Platforms (`supported_formats`)

`tidyproteomics` provides built-in parsers for major quantitative proteomics software packages, as well as an extensible regular-expression engine for custom flat tables (`.csv`, `.tsv`, `.xlsx`).

### Platform Support Matrix

| Platform | Analyte Level | Expected Input Files | Config File Name | Notes |
| :--- | :--- | :--- | :--- | :--- |
| **ProteomeDiscoverer** | `proteins`, `peptides` | `*.xlsx` (Level 1 export) | `ProteomeDiscoverer_proteins.tsv`, `ProteomeDiscoverer_peptides.tsv` | Requires raw abundance columns (e.g., `Abundance: F1: Sample, Control`); extracts sample/file metadata via regex. |
| **MaxQuant** | `proteins` | `proteinGroups.txt` | `MaxQuant_proteins.tsv` | Imports `Intensity` or `LFQ intensity` columns; parses protein IDs and PSM counts. |
| **MaxQuant** | `peptides` | `evidence.txt` | `MaxQuant_peptides.tsv` | Imports peptide-level evidence and intensities across runs. |
| **FragPipe** | `proteins` | `combined_protein.tsv` | `FragPipe_proteins.tsv` | Reads MaxLFQ intensities, spectral counts, and annotations. |
| **FragPipe** | `peptides` | `combined_peptide.tsv` | `FragPipe_peptides.tsv` | Reads peptide-level intensities across experiments. |
| **DIA-NN** | `peptides` | `*.tsv` (report.tsv) | `DIA-NN_peptides.tsv` | Reads `Ms1.Area` or precursor/peptide quantities across runs. |
| **Skyline** | `peptides` | `*.csv` (MSstats peptide report) | `Skyline_peptides.tsv` | Requires standard Skyline MSstats export format. |
| **mzTab** | `proteins`, `peptides` | `*.mzTab` (v1.0.0) | `mzTab_proteins.tsv`, `mzTab_peptides.tsv` | Standard HUPO-PSI tabular format. |
| **Custom / User-Defined** | `proteins`, `peptides` | `*.csv`, `*.tsv`, `*.xlsx` | Custom TSV configuration table | Fully customizable via column mappings and regex extraction. |

### Configuration File Architecture
The import mechanism maps tabular data into the `tidyproteomics` schema using a configuration TSV table:
- `category`: Must include `sample`, `identifier`, `quantitative`, and optionally `accounting`, `annotation`, `impute`, `filter`.
- `column_defined`: Internal destination field (e.g., `sample`, `sample_file`, `protein`, `peptide`, `abundance_raw`, `match_between_runs`, `num_peptides`).
- `column_import`: Exact or regex matching string for the column header in the raw file.
- `pattern_extract`: Regex with lookaheads/lookbehinds to extract values from headers (e.g., extracting sample names from `Abundance: F1: Sample, WT`).
- `pattern_remove`: Regex pattern matching substrings to strip from values.
- `pattern_split`: Character or regex delimiter to split multiple IDs (e.g., `;` for protein accessions).
- `pivot`: `TRUE` if quantitative data is in wide format across samples and requires melting into long format.

### Exporting and Customizing Configs
```r
# Export a default config to the local working directory to inspect or edit
export_config(platform = "ProteomeDiscoverer", analyte = "proteins")
export_config(platform = "MaxQuant", analyte = "peptides")

# Import data using a modified custom configuration table
my_data <- import(
  files = "custom_data_export.tsv",
  platform = "CustomPlatform",
  analyte = "proteins",
  path = "./custom_config.tsv"
)
```

### Match-Between-Runs (MBR)
When MBR is enabled in upstream search software (e.g., MaxQuant, ProteomeDiscoverer), imputed or transferred identifications can be flagged or removed:
```r
# Remove match-between-runs observations
data_no_mbr <- rm.mbr(my_data)
```

---

## 3. Complete API Reference (`api_reference`)

### 3.1 Data Import, Export & Local Persistence

#### `import(files, platform, analyte, path = NULL)`
Imports raw output files into a `tidyproteomics` S3 object.
- **`files`**: Character vector of file paths to import.
- **`platform`**: Name of the platform (`"ProteomeDiscoverer"`, `"MaxQuant"`, `"FragPipe"`, `"DIA-NN"`, `"Skyline"`, `"mzTab"`).
- **`analyte`**: `"proteins"` or `"peptides"`.
- **`path`**: (Optional) Path to a user-defined TSV configuration file.
- **Returns**: A `tidyproteomics` list object.

#### `export_quant(data, file_name = NULL, raw_data = TRUE, normalized = FALSE, scaled = c('none', 'between', 'proportion'))`
Exports the primary quantitative abundance matrix as a tibble or table file (`.xlsx`, `.csv`, `.tsv`, `.rds`).
- **`data`**: `tidyproteomics` data object.
- **`file_name`**: Optional file path. If `NULL`, returns a tibble.
- **`raw_data`**: Logical. Include raw abundances (default `TRUE`).
- **`normalized`**: Character string specifying normalization method to export (or `FALSE`).
- **`scaled`**: Scaling method: `'none'`, `'between'`, or `'proportion'`.
- **Returns**: Tibble (if `file_name` is NULL) or writes file and returns tibble.

#### `export_analysis(data, ..., .analysis, file_name = NULL)`
Exports statistical analysis results (differential expression, enrichment).
- **`data`**: `tidyproteomics` data object.
- **`...`**: Two-sample comparison expression (e.g., `knockdown/control`).
- **`.analysis`**: `"expression"` or `"enrichment"`.
- **`file_name`**: Optional destination file name (`.xlsx`, `.csv`, `.tsv`, `.rds`).
- **Returns**: Tibble of analysis results.

#### `export_compexp(data, ..., file_name = NULL)`
Exports comparative expression data between multiple pairwise comparisons.

#### `export_config(platform, analyte = c('proteins', 'peptides'))`
Copies the internal configuration TSV file for a given platform and analyte to the current working directory.

#### `export_llm_skill(destination = ".", file_name = "SKILL.md", overwrite = TRUE)`
Exports this comprehensive `SKILL.md` reference guide to the local working directory or specified path.
- **`destination`**: Destination directory (default `"."`).
- **`file_name`**: Output filename (default `"SKILL.md"`).
- **`overwrite`**: Logical; overwrite existing file (default `TRUE`).
- **Returns**: Invisible path to the saved file.

#### `prompt_template(type = c("discovery", "profiling", "pairwise", "multicontrast", "extracting"))`
Returns standardized, high-quality prompt templates designed to accompany `SKILL.md` for LLM-driven post-analysis workflows (see Section 6).
- **`type`**: Template flavor to generate:
  - `"discovery"`: Comprehensive Phase 2 biological discovery prompt (*"What is interesting about this dataset?"*).
  - `"profiling"`: Phase 1 reconnaissance prompt generating a minimal R script to extract sample names, counts, and annotations.
  - `"pairwise"`: Targeted two-sample differential abundance prompt.
  - `"multicontrast"`: Batch contrast matrix and multi-term pathway enrichment prompt.
  - `"extracting"`: Venn/Euler diagram overlap and segment extraction prompt.
- **Returns**: Prints template to console and invisibly returns the character string.

#### `load_local(file_name)` / `save_local(data, file_name)`
Loads or saves a `tidyproteomics` object locally as an RDS file with tracking.

#### `save_table(data, file_name)`
Utility function to export any data frame or tibble to `.csv`, `.tsv`, or `.xlsx`.

#### `path_to_package_data(item)`
Returns the absolute path to bundled package example datasets (`"proteins"`, `"peptides"`, `"fasta"`, `"SKILL"`).

#### `extract(data, values = NULL, na.rm = FALSE)`
Extracts quantitative abundance data from the internal `$quantitative` store into a clean, tidy long-format `tibble`. Designed for direct interoperability with standard `tidyverse` pipelines (`dplyr`, `ggplot2`, `readr::write_csv`).
- **`data`**: `tidyproteomics` data object.
- **`values`**: Character vector of quantitative columns to extract (e.g., `"raw"`, `"median"`, `"linear"`, `"loess"`). If `NULL` (default), extracts the active quantitative source (`data$quantitative_source`).
  - *Single method*: Returns columns `sample`, `replicate`, identifier (`protein` or `peptide`), `abundance`, and `origin`.
  - *Multiple methods* (e.g., `values = c("raw", "median", "loess")`): Melts all requested normalizations into a single long tibble with an ordered factor column `origin`, ideal for direct visual comparison of normalization effects.
- **`na.rm`**: Logical (default `FALSE`). When `TRUE`, removes missing or non-positive abundance values (`!is.na(abundance) & abundance > 0`).
- **Returns**: A tidy `tibble` of quantitative observations.
```r
# Extract active quantitative abundance
tbl_quant <- data %>% extract()

# Extract multiple normalizations for comparative visualization
tbl_comp <- data %>% extract(values = c("raw", "median", "loess"), na.rm = TRUE)

# Extract and export to CSV
data %>% extract(na.rm = TRUE) %>% readr::write_csv("tidy_abundance.csv")
```

#### `as.data.frame(x, row.names = NULL, optional = FALSE, shape = c('long', 'wide'), values = NULL, drop = NULL, ...)` (S3 method)
Converts a `tidyproteomics` data object into a standard `tibble` / `data.frame` for seamless integration into custom R workflows.
- **`shape`**: Format of the quantitative matrix:
  - `'long'` (default): One observation per row with sample metadata, accounting, and annotations joined.
  - `'wide'`: One row per protein/peptide with sample quantities spread across columns.
- **`values`**: Normalization method to include (e.g. `'raw'`, `'median'`, `'loess'`). Defaults to the active quantitative source (`data$quantitative_source`).
- **`drop`**: Character vector of accounting or annotation columns to exclude from the resulting data frame.
```r
# Convert to long tibble
df_long <- data %>% as.data.frame() %>% tibble::as_tibble()

# Convert to wide table dropping selected annotation columns
df_wide <- data %>% as.data.frame(shape = "wide", drop = c("description", "biological_process"))
```


---

### 3.2 Data Inspection, Summaries & Helper Utilities

#### `print(x, ...)` (or entering `data` into console)
Exposes a concise, high-level structural summary of the `tidyproteomics` data object. Enforced as the first inspection step upon importing or transforming data:
- **Data Contents**: Origin software (e.g. MaxQuant, ProteomeDiscoverer), analyte level (`proteins` or `peptides`), file count, memory footprint, log10 dynamic range, and % missing values.
- **Sample Names**: Total count and a comma-separated list of all unique sample names.
- **Number of Proteins / Peptides**: Total unique quantified identifiers (e.g., `2450 proteins`).
- **Imputed Method**: If missing values were imputed, explicitly indicates the method used (e.g., `*imputed via randomforest`, `*imputed via min`).
- **Automatic Normalization Selection**: Displays the active quantitative source selected automatically by `select_normalization()` or manually set (e.g., `*normalized loess`).
- **Annotations**: Summary of all biological annotation categories present in the object (e.g., `(4) biological_process cellular_component molecular_function reactome_pathway`).
- **Analyses**: Any completed differential expression comparisons and pathway enrichments.

#### `summary(object, by = c("global", "sample", "normalization", ...), contamination = NULL, ...)`
Primary function for comprehensive **data quality checks** and **normalization checks**. Computes summary statistics returned as a tibble (or printed to console) as a passthrough function without mutating the underlying object:
- **Data Quality Check (`by = "sample"` or `by = "global"`)**:
  - `summary(data)`: Global data quality check reporting total proteins, peptides, peptides per protein, quantifiable fraction, median CV, and PCA dimensionality (PC1 and PC2 variance explained).
  - `summary(data, "sample")`: Per-sample data quality check reporting sample-level identifications, quantifiable fraction, replicate CV, and PCA coordinates (PC1, PC2) to spot outlier runs or sample loading inconsistencies.
  - `summary(data, contamination = 'CRAP')`: Contaminant data quality check screening for common mass spectrometry contaminants against a regular expression pattern (e.g., `'CRAP'`, `'Keratin'`, `'Trypsin'`, `'BSA'`).
    > [!IMPORTANT]
    > **Always Use Named Argument `contamination = 'CRAP'`**:
    > To summarize contamination, always call `summary(data, contamination = 'CRAP')`. **Never write `summary(data, "contamination")`**! Passing `"contamination"` positionally sets `by = "contamination"`, which causes `tidyproteomics` to look for a non-existent column named `"contamination"` and returns `NULL`.
- **Normalization Check (`by = "normalization"`)**:
  - `summary(data, "normalization")`: Essential normalization check comparing all candidate normalized datasets (raw, median, linear, loess, randomforest, etc.) across identification recovery, quantifiable fraction, variance reduction (median CV), and PC1 & PC2 variance explanation, allowing users and agents to rigorously validate automatic normalization selection before downstream differential analysis.


#### `operations(data)`
Prints and returns the full sequential history of operations performed on the dataset.

#### Data Object Inspection Getters (`get_*()`)
Exposed inspection helpers to examine variables, levels, metrics, and annotations within a `tidyproteomics` data object:

- **`get_sample_names(data)`**: Returns a character vector of all unique sample condition names.
- **`get_quant_names(data)`**: Returns a character vector of available quantitative data columns (e.g., `"raw"`, `"median"`, `"linear"`, `"loess"`).
- **`get_accountings(data)`**: Returns a character vector of accounting column names (e.g., `"num_peptides"`, `"num_unique_peptides"`).
- **`get_annotation_terms(data)`**: Returns a character vector of available annotation category terms (including the primary identifier, e.g., `"protein"`, `"gene_name"`, `"description"`).
- **`get_annotations(data, term = NULL)`**: Returns a 2-column tibble (`identifier`, `annotation`) for the specified annotation term.
- **`show_annotations(data)`**: Displays a formatted console summary of available annotation categories.
- **`get_variables(data, segment = c("experiments", "quantitative", "annotations", "accounting"))`**: Returns a character vector of variable names within the requested object segment.
- **`get_unique_variables(data, variable = NULL)`**: Returns a vector of unique values for any specified variable across the object (including expanded annotation values).


#### `annotate(data, annotations = NULL, duplicates = c("replace", "merge", "leave"))`
Attaches or updates biological annotations (e.g., Gene Ontology terms, Reactome pathways, updated gene symbols/descriptions) in a `tidyproteomics` data object.
- **`data`**: `tidyproteomics` data object.
- **`annotations`**: Data frame / tibble containing the identifier column (`protein` or `peptide`), `term` (character string naming the annotation category, e.g. `'biological_process'`, `'molecular_function'`, `'gene_name'`), and `annotation` (character string describing the term/value).
- **`duplicates`**: How to handle duplicate terms for the same identifier:
  - `"replace"` (default): Keeps the newest annotation value, replacing existing ones.
  - `"merge"`: Concatenates multiple annotation strings using `"; "` delimiters.
  - `"leave"`: Retains multiple rows without merging or discarding.
- **Returns**: Updated `tidyproteomics` data object.
```r
# Add custom GO annotations formatted as identifier, term, and annotation
data <- data %>% annotate(custom_go_table, duplicates = "merge")
```

#### `experimental(data)` / `experimental_groups(data)`
Inspects or assigns experimental factors and groupings.

#### `reassign(data, ..., .replace = NULL)`
Modifies sample labels, groupings, or metadata mappings.

#### `origins(data, ..., .replace = NULL)`
Assigns or modifies biological replicate or subject origin (`sample_origin`) in `data$experiments` for each input file/run. Used critically by `expression(..., .method = "mixed-effects")` / `expression_me()` to account for repeated measures, multi-region sampling, or subject batches in linear mixed-effects models:
```r
data <- data %>%
  origins(sample_file == "F1", .replace = "donor_1") %>%
  origins(sample_file == "F2", .replace = "donor_1") %>%
  origins(sample_file == "F3", .replace = "donor_2")
```

#### `merge(data1, data2, ...)` / `merge_tidyproteomics(data_list, ...)`
Combines multiple `tidyproteomics` data objects into a single unified object via S3 `merge()` dispatch or `merge_tidyproteomics()`.

#### `rm.mbr(data)`
Removes quantitative values flagged as Match-Between-Runs (MBR).

#### `hdf(x, n = 5)`
Convenience helper returning `head(x, n)` converted to a standard `data.frame` for clean printing.

#### `%like%` (Regex Matching Operator)
Infix operator for case-insensitive regex matching. Syntax: `string %like% "pattern"`. Widely used inside `subset()`:
```r
data %>% subset(description %like% "Ribosome")
data %>% subset(!description %like% "Keratin")
```

---

### 3.3 Subsetting & Filtering

#### `subset(data, ..., .verbose = TRUE)`
Filters proteins or peptides based on annotations, abundance thresholds, or identification metrics.
```r
# Filter by annotation regex
data_ribo <- data %>% subset(description %like% "Ribosome")

# Remove contaminants
data_clean <- data %>% subset(!description %like% "Keratin|Trypsin|Albumin")

# Filter by accounting (e.g., minimum unique peptides)
data_high_conf <- data %>% subset(num_unique_peptides >= 2)
```

#### `intersection(data, .include = NULL, .exclude = NULL)`
Extracts the underlying data corresponding to specific Venn and Euler diagram segments (`extracting.Rmd`). While `plot_venn()` and `plot_euler()` visually highlight overlapping and distinct proteins between sample groups, `intersection()` extracts the actual data subset belonging to any given diagram segment.
- **`.include`**: Character vector of sample groups that *must* contain the protein.
- **`.exclude`**: Character vector of sample groups that *must not* contain the protein.
- **Key S3 Object Architecture Feature**:
  `intersection()` returns a **valid, fully synchronized `tidyproteomics` data object**, not just a vector of IDs. It automatically subsets `$quantitative`, `$accounting`, `$annotations`, and `$experiments` simultaneously, updating the object so that only the samples and raw files associated with the extracted segment remain (e.g., removing unused sample files).
- **Core Use Cases**:
  - *Condition-Exclusive Biomarkers*: Isolate proteins uniquely identified in one group (e.g. induced in treatment or lost upon disease/knockout).
  - *Complex Multi-Group Overlaps*: Extract specific intersections across 3, 4, or more conditions (e.g., proteins present in groups A and B, but absent in C and D).
  - *Pipeable Downstream*: Because it returns a `tidyproteomics` object, you can immediately inspect it with `print()`, check it with `summary()`, or extract a clean tabular data frame using `extract(na.rm = TRUE)`.

```r
# 1. Condition-exclusive extraction (proteins unique to control, absent in knockdown):
control_exclusive <- hela_proteins %>%
  subset(imputed == 0) %>%
  intersection(.include = c('control'), .exclude = c('knockdown'))

# 2. Shared core proteome (proteins present in both conditions):
core_shared <- hela_proteins %>%
  subset(imputed == 0) %>%
  intersection(.include = c('control', 'knockdown'))

# 3. Complex 4-way Venn segment (proteins in a, c, and d, but NOT in b):
exp_proteins %>%
  intersection(.include = c('a', 'c', 'd'), .exclude = c('b'))

# 4. Pull underlying tidy table of exclusive proteins:
tbl_exclusive <- control_exclusive %>% extract(na.rm = TRUE)
```


---

### 3.4 Exploratory Data Analysis & Visualization

All `plot_*()` functions in `tidyproteomics` integrate directly with the internal saving engine `plot_save()`. They do **not** require external plotting libraries or `ggplot2::ggsave()`.

#### Built-in Saving Inherited from `plot_save()`
Every `plot_*()` function inherits the following parameters from `plot_save()` via `...`:
- **`destination`**: character string specifying output behavior:
  - `"plot"` (default): Returns the ggplot2 / grob object directly to the R console or plotting window for interactive inspection.
  - `"save"` or `"png"`: Saves the plot to disk as a PNG file.
  - `"svg"`: Saves as an SVG vector file.
  - `"tiff"`: Saves as a high-resolution TIFF file.
  - `"jpeg"`: Saves as a JPEG file.
- **`height`**: Numeric plot height in inches (default `5`, allowed range `1` to `30`).
- **`width`**: Numeric plot width in inches (default `8`, allowed range `1` to `30`).
- **`...`**: Additional parameters forwarded to `ggplot2::ggsave()` (e.g. `dpi = 300`).

Output filenames are constructed automatically using the pattern: `plot_{file_name}.{destination}` (e.g., `plot_proteins_counts.png`, `plot_proteins_pca.png`, `plot_proteins_volcano_knockdown_control.png`).

#### Pipeline Chaining & Console Printing Behavior
When `destination` is set to save (`"save"`, `"png"`, `"svg"`, `"tiff"`, `"jpeg"`), `plot_*()` writes the graphic to disk and **returns the original `tidyproteomics` data object**.

This design allows you to concatenate multiple plot commands sequentially in a clean, readable pipeline:
```r
# Concatenate multiple QC plots in a single pipeline
data %>%
  plot_counts(destination = "png") %>%
  plot_variation_cv(destination = "png") %>%
  plot_dynamic_range(destination = "png") %>%
  plot_pca(destination = "png") -> data_qc
```

> [!WARNING]
> **Capture the Return Value to Avoid Unintended Printing**:
> Because saving returns the `tidyproteomics` data object, running an unassigned pipeline interactively causes R to evaluate and print the returned object, which invokes `print(data)` and displays the full quantitative proteomics summary in the console.
> Always **assign the result** (e.g., `data <- data %>% plot_...`) or wrap in **`invisible()`**:
> ```r
> invisible(
>   data %>%
>     plot_counts(destination = "png") %>%
>     plot_pca(destination = "png")
> )
> ```

> [!CAUTION]
> **STRICT RULE FOR LLMs: NEVER USE `ggsave()` WITH `plot_*()`**:
> LLMs frequently make the mistake of piping `plot_*()` into `ggplot2::ggsave()` or writing `ggsave(filename, data %>% plot_*())`.
> **This will fail with a fatal runtime error!**
> When saving is enabled, `plot_*()` returns a `tidyproteomics` S3 object, which `ggsave()` immediately rejects because it expects a ggplot object, not a `tidyproteomics` object.
> Saving is handled entirely internally by specifying `destination = "png"` (or `"save"`, `"svg"`, etc.) directly inside `plot_*()`.

#### Available Plotting Functions
- **`plot_counts(data, show_replicates = TRUE, impute_max = 0.5, palette = 'YlGnBu', ...)`**: Bar chart of observed proteins/peptides per sample and replicate.
- **`plot_quantrank(data, palette = 'BuGn', ...)`**: Rank-abundance curve showing dynamic range and abundance distribution.
- **`plot_variation_cv(data, palette = 'YlGnBu', ...)`**: Distribution of Coefficients of Variation (%CV) across replicates.
- **`plot_variation_pca(data, ...)`**: Scree / variance explained plot for Principal Component Analysis.
- **`plot_pca(data, palette = 'Set1', ...)`**: 2D PCA score plot colored by sample groups with confidence ellipses.
- **`plot_dynamic_range(data, palette = 'YlGnBu', ...)`**: Boxplot/violin plot of dynamic range across samples.
- **`plot_venn(data, ...)`** / **`plot_euler(data, ...)`**: Venn and Euler diagrams showing overlapping identifications between sample groups (`extracting.Rmd`). Commonly paired with `subset(imputed == 0)` to inspect true observed identifications, followed by `intersection()` to extract the underlying proteins belonging to specific diagram regions.
- **`plot_heatmap(data, ...)`**: Clustered heatmap of protein abundances across samples.
- **`plot_proportion(data, ...)`**: Stacked proportion visualization of abundance across samples.
- **`plot_save(plot, data, file_name, destination = c("plot", "save", "png", "svg", "tiff", "jpeg"), height = 5, width = 8, ...)`**: The core underlying saving engine.


---

### 3.5 Normalization

#### `normalize(data, ..., .method = c("scaled", "median", "linear", "limma", "loess", "svm", "randomforest"), .cores = 1)`
Applies one or more normalization algorithms to the raw quantitative data.
- **Methods**:
  - `"median"`: Global median alignment across samples.
  - `"scaled"`: Proportional total signal scaling.
  - `"linear"`: Linear regression-based normalization across intensity range.
  - `"loess"`: Local regression normalization (handles non-linear systematic biases).
  - `"limma"`: Cyclic loess / quantile-based normalization via the `limma` package.
  - `"svm"`: Support Vector Machine regression normalization.
  - `"randomforest"`: Random Forest regression normalization.
  - `"tmt-bridge"`: Normalization using a common reference bridge channel in TMT workflows.
- **Targeted Normalization via Subset**: Normalizes based on a stable subset (e.g., housekeeping proteins or spike-ins) and applies the transformation globally:
```r
data %>% normalize(description %like% "Histone", .method = "median")
```

#### `select_normalization(data, normalization = NULL)`
Evaluates all computed normalization methods based on three quantitative metrics:
1. Minimizing intra-condition CVs.
2. Minimizing PC1 variance driven by technical artifacts.
3. Maximizing preserved dynamic range.
If `normalization` is specified (e.g., `"median"`), manually sets `data$quantitative_source`.

#### `plot_normalization(data)`
Generates comparison plots evaluating raw vs. normalized data distributions.

---

### 3.6 Missing Value Imputation

#### `impute(data, .function = base::min, method = c("row", "column", "matrix"), group_by_sample = FALSE, cores = 2)`
Imputes missing quantitative observations.
- **`.function`**: Imputation function (e.g., `base::min`, `stats::median`, `stats::mean`, `impute.randomforest`).
- **`method`**:
  - `"row"`: Imputes using values of the same feature across other samples.
  - `"column"`: Imputes using values within the same sample distribution (e.g., left-censored minimum).
  - `"matrix"`: Global matrix imputation (e.g., Random Forest across rows and columns).
- **`group_by_sample`**: Logical. Group within sample replicates to prevent cross-condition leakage.

---

### 3.7 Two-Sample Differential Expression Analysis

#### `expression(data, ..., .method = "limma", .p.adjust = "BH", .pairs = NULL)`
Computes differential protein/peptide abundance between two experimental groups.
- **`...`**: Two-sample comparison in ratio syntax: `treatment/control` or `knockdown/wildtype`.
- **`.method`**: Statistical test:
  - `"limma"`: Empirical Bayes moderation via Limma (recommended for proteomics).
  - `"mixed-effects"`: Linear mixed-effects model via `MSstats` (`expression_me()`). Uses `sample_origin` from `data$experiments` (set via `origins()`) as the biological subject mapping (`BioReplicate`).
  - `stats::t.test`: Standard Student's or Welch's t-test.
  - `stats::wilcox.test`: Non-parametric Wilcoxon / Mann-Whitney test.
  - Custom functions: Any two-sample function taking `(x, y)` and returning a p-value.
- **`.p.adjust`**: Correction for multiple testing (default `"BH"` / Benjamini-Hochberg; also `"bonferroni"`, `"holm"`, `"fdr"`).
- **`.pairs`**: List of 2-element character vectors for multiple pairwise comparisons:
```r
data %>% expression(.pairs = list(c("knockdown", "control"), c("drug", "control")))
```

#### `analyze_expressions(data, labels_column = NULL, log2fc_min = 1, significance_max = 0.05, ...)`
Automates end-to-end report generation for all expression analyses stored in `data$analysis` (`automating.Rmd`). Iterates through all pairwise comparisons (previously computed via `expression(treatment/control)` or `expression(.pairs = pairs)`), generating publication-ready volcano plots and tabular exports directly to disk.
- **`data`**: `tidyproteomics` data object with completed `expression()` computations.
- **`labels_column`**: Character name of the annotation column used to label top significant hits (e.g., `'gene_name'` or `'protein'`).
- **`log2fc_min`**: Minimum absolute log2 fold-change threshold for highlighting (default `1`).
- **`significance_max`**: Maximum significance cutoff for highlighting (default `0.05`).
- **`significance_column`**: Significance column to filter (`'adj_p_value'` [default] or `'p_value'`).
- **`show_pannels`**: Logical (default `TRUE`). Displays colored background panels distinguishing up- and down-regulated regions.
- **`show_lines`**: Logical (default `TRUE`). Draws threshold dashed lines for fold-change and significance.
- **`show_fc_scale`**: Logical (default `TRUE`). Displays a secondary linear fold-change axis for biological interpretability.
- **Auto-Generated Artifacts**:
  - `table_proteins_expression_<contrast>.csv`: Full tabular differential abundance matrix for each comparison.
  - `plot_volcano.png` (or contrast-specific): Styled volcano plots saved directly to the local directory.
- **Returns**: Invisibly returns the input `data` object, enabling seamless pipeline chaining.

#### `plot_volcano(data, ..., log2fc_min = 1, significance_max = 0.05, labels = TRUE)`
Produces publication-quality volcano plots (`log2(Fold Change)` vs. `-log10(p-value)` or `q-value`) with automatic ggrepel text labels.
```r
data %>% plot_volcano(knockdown/control, log2fc_min = 1.5, significance_max = 0.01)
```

#### `plot_expression(data, ..., log2fc_min = 1, significance_max = 0.05, top_n = 10, regulation = c("both", "up", "down"), labels_column = "protein")`
Produces faceted boxplots with overlaid individual replicate observations (filled circles for observed, hollow circles for imputed) displaying log2 abundance distributions for the top $N$ differentially expressed proteins between two conditions. Each facet panel includes an indicator crossbar and statistical significance annotation (`*`, `**`, and p-value).
- **`...`**: Two-sample contrast formula (e.g. `knockdown/control`).
- **`log2fc_min` / `significance_max`**: Fold-change and significance cutoff thresholds for candidate filtering.
- **`top_n`**: Number of top candidates to plot post-filtering (default `10`).
- **`regulation`**: Direction of change to include: `"both"` (default), `"up"`, or `"down"`.
- **`labels_column`**: Column name used for facet headers (default `"protein"`, or any annotation term such as `"gene_name"`).
- **`destination`**: Output format (`"plot"`, `"save"`, `"png"`, `"svg"`, `"tiff"`, `"jpeg"`), adhering to the `plot_save()` framework.

```r
data %>%  plot_expression(knockdown/control,
                  log2fc_min = 0.5,
                  significance_column = "p_value",
                  top_n = 6,
                  labels_column = "gene_name",
                  destination = "png")
```

#### `plot_compexp(data, ...)`
Comparative scatter plot contrasting fold changes between two distinct comparisons (e.g., Treatment A vs. Treatment B against common Control).

---

### 3.8 Term & Pathway Enrichment Analysis

#### `enrichment(data, ..., .method = c("gsea", "wilcoxon", "fishersexact"), .terms = "biological_process", .cpu_cores = 1, .pairs = NULL)`
Tests for over-representation or rank enrichment of functional annotation terms.
- **`.method`**:
  - `"gsea"`: Gene Set Enrichment Analysis using fast `fgsea` (ranked by fold change / signed statistic).
  - `"wilcoxon"`: Wilcoxon rank-sum test on protein fold changes.
  - `"fishersexact"`: 2x2 contingency table test on significant vs. non-significant proteins.
- **`.terms`**: Single term or character vector of multiple annotation categories:
  ```r
  .terms = c("biological_process", "cellular_component", "molecular_function", "wiki_pathway", "reactome_pathway")
  ```
- **`.pairs`**: List of 2-element character vectors for automated batch testing across multiple comparisons.

#### `analyze_enrichments(data, significance_max = 0.05, top_n = 50, enriched_up_color = 'blue', enriched_down_color = 'red', ...)`
Automates end-to-end pathway report generation across all pairwise comparisons and annotation terms stored in `data$analysis` (`automating.Rmd`). Iterates through all computed enrichments, outputting CSV data tables and a consolidated summary figure.
- **`data`**: `tidyproteomics` data object with completed `enrichment()` calculations.
- **`significance_max`**: Adjusted p-value cutoff for significance highlighting (default `0.05`).
- **`top_n`**: Number of top terms to display per contrast in the summary plot (default `50`).
- **`enriched_up_color`** / **`enriched_down_color`**: Colors assigned to positive (up) and negative (down) pathway enrichments.
- **Auto-Generated Artifacts**:
  - `table_proteins_enrichment_<contrast>.csv`: Full enrichment results per contrast and term.
  - `proteins_enrichment_<contrast>.png`: A single composite, multi-term horizontal bar plot summarizing all tested annotation categories side-by-side, color-coded by up/down regulation and significance.
- **Returns**: Invisibly returns the input `data` object, allowing it to terminate an automated pipeline.


#### `plot_enrichment(data, ...)`
Generates bubble plots or bar charts displaying normalized enrichment scores (NES) and adjusted p-values.

---

### 3.9 Peptide-to-Protein Collapsing & Protein Sequence Mapping

#### `collapse(data, collapse_to = "protein", assign_by = c("all-possible", "razor-local", "razor-global", "non-homologous"), top_n = Inf, fasta_path = NULL, split_abundance = FALSE, .function = fsum)`
Aggregates peptide-level quantitative data into protein-level objects.
- **`collapse_to`**: Annotation target (default `"protein"`; can also be `"gene_name"`).
- **`assign_by`**:
  - `"razor-local"`: Assigns shared peptides to the protein of highest likelihood within a sample group.
  - `"razor-global"`: Assigns shared peptides based on overall highest abundance protein globally.
  - `"non-homologous"`: Restricts aggregation strictly to unique, non-shared peptides.
  - `"all-possible"`: Contributes peptide abundance to all matching proteins.
- **`.function`**: Fast aggregation functions: `fsum` (default), `fmedian`, `fmean`, `fgeomean`, `fmin`.

#### `protein_map(data, fasta_path)`
Maps peptide sequences and post-translational modifications (PTMs) back to full-length protein sequences parsed from a FASTA file.

#### `plot_protein(mapped_data, protein, ...)`
Visualizes protein sequence coverage, peptide alignments, and modification frequencies along the primary sequence.

---

## 4. Canonical End-to-End Workflows

### Workflow 1: Standard Discovery Proteomics Pipeline
```r
library(dplyr)
library(tidyproteomics)

# 1. Import raw data
raw_data <- import(
  files = "proteinGroups.txt",
  platform = "MaxQuant",
  analyte = "proteins"
)

# 2. Inspect dataset & Data Quality Check
print(raw_data)               # Concise summary: contents, samples, protein count, annotations
raw_data %>% summary("sample") # Data quality check: per-sample counts, CVs, and PCA separation
raw_data %>% summary(contamination = "CRAP") # Contaminant quality check
raw_data %>% plot_counts()

# 3. Filter contaminants
clean_data <- raw_data %>%
  subset(!description %like% "CON__|Keratin|Trypsin") %>%
  subset(num_unique_peptides >= 2)

# 4. Normalize, check, and select best method
norm_data <- clean_data %>%
  normalize(.method = c("scaled", "median", "linear", "loess", "limma")) %>%
  select_normalization()

norm_data %>% summary("normalization") # Normalization check: compare CV reduction and PC separation
print(norm_data)                       # Verifies automatic normalization selection (e.g. *normalized loess)

# 5. Impute missing values (left-censored min within sample group)
imputed_data <- norm_data %>%
  impute(.function = base::min, method = "column", group_by_sample = TRUE)

print(imputed_data)                    # Verifies imputed method (e.g. *imputed min)

# 6. Two-sample differential expression
diff_data <- imputed_data %>%
  expression(knockdown/control, .method = "limma")

# 7. Visualize results (chain multiple plots with built-in saving; capture to prevent printing data object)
diff_data <- diff_data %>%
  plot_volcano(knockdown/control, log2fc_min = 1, significance_max = 0.05, destination = "png") %>%
  plot_pca(destination = "png")

# 8. Export quantitative matrix & analysis tables
diff_data %>% export_quant(file_name = "normalized_abundance.xlsx", normalized = "loess")
diff_data %>% export_analysis(knockdown/control, .analysis = "expression", file_name = "de_results.xlsx")

# 9. Review audit trail
diff_data %>% operations()
```

### Workflow 2: Peptide-to-Protein Collapsing Pipeline
```r
library(dplyr)
library(tidyproteomics)

# 1. Import peptide-level data
peptides <- import(
  files = "evidence.txt",
  platform = "MaxQuant",
  analyte = "peptides"
)

# 2. Collapse peptides to proteins using razor-global logic
proteins <- peptides %>%
  collapse(
    collapse_to = "protein",
    assign_by = "razor-global",
    fasta_path = "uniprot_human.fasta",
    .function = fsum
  )

# 3. Proceed with protein-level analysis (save volcano directly via destination = 'png')
results <- proteins %>%
  normalize(.method = "median") %>%
  expression(treatment/control) %>%
  plot_volcano(treatment/control, destination = "png")
```

### Workflow 3: Batch Automation Pipeline (`automating.Rmd`)
```r
library(dplyr)
library(tidyproteomics)

# 1. Define contrast pairs and functional annotation terms
pairs <- list(
  c("knockdown", "control"),
  c("rescue", "control"),
  c("rescue", "knockdown")
)

terms <- c(
  "biological_process",
  "cellular_component",
  "molecular_function",
  "wiki_pathway",
  "reactome_pathway"
)

# 2. Execute end-to-end automated analysis in a single fluent pipeline:
#    - expression(.pairs = pairs): batch differential testing for all pairs
#    - enrichment(.pairs = pairs, .terms = terms): multi-term GSEA for all pairs
#    - analyze_expressions(): auto-generates styled volcano plots & CSV tables for all pairs
#    - analyze_enrichments(): auto-generates multi-term horizontal bar plots & CSV tables
results <- my_data %>%
  expression(.pairs = pairs, .method = "limma") %>%
  enrichment(.pairs = pairs, .terms = terms, .method = "wilcoxon", .cpu_cores = 4) %>%
  analyze_expressions(labels_column = "gene_name", log2fc_min = 1, significance_max = 0.05) %>%
  analyze_enrichments(significance_max = 0.05, top_n = 50)

# 3. Auto-generated files written directly to disk:
# - table_proteins_expression_<contrast>.csv (differential expression tables)
# - plot_volcano.png (styled volcano plots with thresholds and secondary fold-change scales)
# - table_proteins_enrichment_<contrast>.csv (pathway enrichment tables across all terms)
# - proteins_enrichment_<contrast>.png (composite multi-term horizontal bar plot summarizing all tested terms)
```

### Workflow 4: Venn Diagramming & Underlying Segment Extraction (`extracting.Rmd`)
```r
library(dplyr)
library(tidyproteomics)

# 1. Visualize proteome overlap across sample groups
# Filter non-imputed proteins (imputed == 0) to evaluate true mass spectrometry identifications
my_data %>%
  subset(imputed == 0) %>%
  plot_venn(destination = "png")

# 2. Extract single-condition exclusive proteins (e.g. Unique to 'control', absent in 'knockdown')
# intersection() trims $quantitative, $accounting, $annotations, and $experiments simultaneously
control_exclusive <- my_data %>%
  subset(imputed == 0) %>%
  intersection(.include = c("control"), .exclude = c("knockdown"))

print(control_exclusive) # Only 1 sample and its associated raw files remain in the object!

# 3. Extract tidy tabular data frame of exclusive proteins with abundances and annotations
tbl_control_exclusive <- control_exclusive %>%
  extract(na.rm = TRUE)

# Export to CSV
readr::write_csv(tbl_control_exclusive, "proteins_exclusive_control.csv")

# 4. Complex multi-group segment extraction (e.g. 4-way comparison)
# Extract proteins present in groups a, c, and d, but absent in b:
acd_proteins <- my_data %>%
  subset(imputed == 0) %>%
  intersection(.include = c("a", "c", "d"), .exclude = c("b"))

# Extract proteins shared between a and b, but absent in c and d, and convert directly to tibble:
ab_shared_tbl <- my_data %>%
  subset(imputed == 0) %>%
  intersection(.include = c("a", "b"), .exclude = c("c", "d")) %>%
  extract(na.rm = TRUE)
```

---

## 5. LLM Prompting & Autonomous Agent Integration Guide

When using an LLM or AI coding agent (e.g., Antigravity, Claude, ChatGPT, Gemini, Copilot) with `tidyproteomics`:

### Guidelines for LLMs
1. **Enforce `print(data)` for Immediate Orientation & State Verification**:
   Always run `print(data)` (or evaluate `data`) immediately after loading a dataset and following transformations (e.g. normalization, imputation). It provides a concise summary of:
   - **Data contents** (origin platform, analyte type, file count, memory size, dynamic range, percent missing values)
   - **Sample names** (essential for specifying contrast formulas without guessing labels)
   - **Number of proteins / identifiers** quantified
   - **Imputed method** (explicitly displays applied imputation algorithm, e.g., `*imputed via randomforest`)
   - **Automatic normalization selection** (verifies active quantitative source, e.g., `*normalized loess`)
   - **Annotations** (identifies loaded biological categories such as `biological_process` or `reactome_pathway`)
2. **Enforce `summary(data)` for Data Quality & Normalization Checks**:
   - **Data Quality Check**: Run `summary(data, "sample")` to inspect replicate counts, missingness, median CV, and PCA separation (PC1, PC2) across samples to catch technical outliers and uneven loading.
   - **Contamination Check**: Run `summary(data, contamination = 'CRAP')` (or `contamination = "Keratin|Trypsin"`) to detect potential sample preparation contaminants before proceeding with differential expression.
     > **CRITICAL**: Always use the named parameter `contamination = 'CRAP'`. Never write `summary(data, "contamination")`, which erroneously searches for a column named `"contamination"` and returns `NULL`.
   - **Normalization Check**: After normalizing, run `summary(data, "normalization")` to evaluate all candidate normalizations against raw data (comparing CV reduction, quantifiable fraction, and PC1 & PC2 variance explanation) to validate automatic normalization selection.
3. **NEVER Call `ggsave()` on `plot_*()` Outputs**:
   - All `plot_*()` functions inherit saving capabilities directly from `plot_save()` via `destination = c("png", "svg", "tiff", "jpeg", "save")`, `height`, and `width`.
   - **Fatal Error Warning**: Calling `ggsave()` on `data %>% plot_*(destination = "png")` fails because saving mode writes the file and returns the `tidyproteomics` data object, not a ggplot object!
   - Always save plots directly by specifying `destination = "png"` (or `"svg"`, `"save"`, etc.) inside the `plot_*()` call itself.
4. **Capture `plot_*()` Pipelines to Prevent Console Printing**:
   - Because `plot_*(..., destination = "png")` returns the `tidyproteomics` data object, multiple plot functions can be concatenated sequentially in a single pipe chain:
     ```r
     data %>%
       plot_counts(destination = "png") %>%
       plot_pca(destination = "png") -> data
     ```
   - **Crucial Console Rule**: If an unassigned pipeline is executed interactively, R automatically evaluates and prints the returned `data` object at the end of the pipe, triggering `print(data)` and dumping the full object summary to the console.
   - Always assign the output (e.g., `data <- data %>% plot_...`) or wrap in `invisible(...)`.
5. **Always Check Object History**: Run `operations(data)` to inspect what normalizations or filters have already been applied before chaining additional modifications.
6. **Use `%like%` for Pattern Matching**: When filtering annotations in `subset()`, use `description %like% "regex"` rather than `grepl()` inside base subset.
7. **Preserve Ratios Syntax**: When defining comparisons in `expression()`, use unquoted ratio notation `numerator/denominator` (e.g. `treatment/control`). The resulting `log2_foldchange` is calculated as $\log_2(\text{numerator}) - \log_2(\text{denominator})$.
8. **Select Normalization Before Statistical Testing**: Always verify that `select_normalization()` has been called or specify which normalization to test, ensuring tests do not run against unnormalized raw abundance unintentionally.
9. **Fast Collapsing**: When aggregating large peptide datasets, prefer `fsum` or `fmedian` over generic R functions, as they are optimized for missing value removal.
10. **Automate Multi-Contrast Screens via `automating.Rmd`**:
   When analyzing multi-condition screens, timepoints, or multi-sample comparisons, NEVER write manual `for` loops or repetitive `lapply` blocks to execute and save differential expressions or enrichments. Instead, pass `.pairs = list(...)` to `expression()`, pass `.pairs` and `.terms = c(...)` to `enrichment()`, and chain `analyze_expressions(labels_column = "gene_name")` and `analyze_enrichments(significance_max = 0.05)`. These automated functions handle iteration internally, produce formatted tables (`table_proteins_expression_*.csv`, `table_proteins_enrichment_*.csv`), render styled volcano plots, and generate consolidated multi-term horizontal bar charts (`proteins_enrichment_*.png`).
11. **Extract Venn/Euler Diagram Segments via `extracting.Rmd`**:
   When users ask "What proteins are unique to group A?", "What proteins are lost upon knockout?", or want underlying data from Venn overlaps:
   - Filter non-imputed proteins with `subset(imputed == 0)` to guarantee true observed identifications.
   - Use `intersection(.include = ..., .exclude = ...)` to isolate the specific diagram segment. Never attempt manual set operations or filtering on data frames.
   - Recognize that `intersection()` returns a valid, fully synchronized `tidyproteomics` S3 object (retaining only the relevant samples and files).
   - Pipe into `extract(na.rm = TRUE)` to produce a clean long-format tibble ready for `readr::write_csv()` or further `dplyr` filtering.
12. **Exporting Skill Locally**: To inspect or store this reference in an agent workspace:
   ```r
   tidyproteomics::export_llm_skill()
   ```

---

## 6. Two-Phase Discovery Strategy & Prompt Templates

When researchers ask an LLM *"Tell me what is interesting about this dataset"*, open-ended queries often produce code with guessed sample names or premature statistical assumptions. The most reliable approach is a **Two-Phase Workflow**:

```
[User has raw data]
       │
       ▼  (Phase 1: Profiling Prompt)
[LLM writes 5-line R reconnaissance script]
       │
       ▼  (User runs script in R, copies summary output back to LLM)
[LLM knows exact sample names, replicate counts, & annotations]
       │
       ▼  (Phase 2: Discovery Prompt: "What is interesting about this dataset?")
[LLM produces complete analysis: QC -> Normalization -> Volcano -> GSEA -> Biological Interpretation]
```

Researchers can generate these prompts directly in R using `tidyproteomics::prompt_template()` or copy the templates below:

### Prompt Template 1: Phase 1 — Dataset Profiling & Reconnaissance
*Use this prompt first to let the LLM generate a quick metadata extraction script:*

```markdown
You are an expert computational proteomics specialist using the R package `tidyproteomics`.
Refer to the instructions in SKILL.md.

I have imported my quantitative mass spectrometry dataset in R, but need to inspect its structure before performing statistical analysis.

Please provide a minimal, executable R script that:
1. Loads `tidyproteomics` and `dplyr`.
2. Imports the data file located at "[PATH_TO_DATA_FILE]" from platform "[PLATFORM: e.g. MaxQuant / ProteomeDiscoverer / FragPipe / DIA-NN]" as "[ANALYTE: proteins / peptides]".
3. Runs `print(data)` to inspect the concise summary of data contents, sample names, number of proteins, imputed method, automatic normalization selection, and annotations.
4. Runs `summary("sample")` as a data quality check to inspect replicate-level counts, CVs, and PCA separation.
5. Prints all unique sample names using `get_sample_names()`.
6. Prints available annotation categories using `show_annotations()`.
7. Checks for potential sample contaminants using `summary(contamination = 'CRAP')`.

I will run this script in R and share the console output with you so we can design the comparative analysis.
```

### Prompt Template 2: Phase 2 — "What is Interesting About This Dataset?" (Discovery Engine)
*Use this prompt once sample names and annotations are known:*

```markdown
You are an expert quantitative proteomics data scientist. You have been provided with the `tidyproteomics` package instructions in SKILL.md.

I have an imported `tidyproteomics` dataset with the following experimental setup:
- Organism: [e.g. Homo sapiens / Mus musculus]
- Biological Question: [e.g. Mechanism of action of Drug X / Knockout phenotype]
- Sample Groups: [e.g. Treatment = "drug_10uM", Control = "dmso", 3 replicates each]
- Available Annotations: [e.g. gene_name, description, biological_process, reactome_pathway]

Your goal is to answer: "What is biologically and technically interesting about this dataset?"

Generate a complete, self-contained R analysis script using `tidyproteomics` that:
1. Evaluates Quality Control & Global Structure (Data Quality Check):
   - Performs a data quality check with `summary("sample")` to inspect replicate CVs, quantifiable fraction, and PCA coordinates.
   - Saves QC plots using built-in `destination = "png"` (e.g. `plot_variation_cv(destination = "png")`, `plot_dynamic_range(destination = "png")`, `plot_pca(destination = "png")`).
   - Chains plot calls directly without `ggsave()`, capturing returned data objects in variables (`data <- data %>% plot_...`) to avoid console printing.
2. Selects Optimal Normalization & Normalization Check:
   - Tests candidate normalizations (`c("median", "linear", "loess", "limma")`) and selects the optimal method (`select_normalization()`).
   - Runs `summary("normalization")` as a normalization check to confirm CV reduction and variance preservation across candidate normalizations.
   - Verifies the automatic normalization selection via `print(data)`.
3. Handles Missing Values:
   - Imputes missing data appropriately (e.g. left-censored minimum within sample groups or random forest).
   - Verifies the imputed method via `print(data)`.
4. Uncovers Differential Expression:
   - Compares Treatment vs Control using `expression(Treatment/Control, .method = "limma")`.
   - Generates and saves a publication-ready volcano plot (`plot_volcano(destination = "png")`, log2fc_min = 1.0, significance_max = 0.05).
5. Explores Biological Mechanism (Pathway Enrichment):
   - Runs GSEA (`enrichment(Treatment/Control, .method = "gsea", .terms = "biological_process")`).
   - Generates and saves an enrichment plot (`plot_enrichment(destination = "png")`).
6. Highlights Key Findings:
   - Extracts and prints the top 10 most significantly upregulated and downregulated proteins.
   - Extracts and prints the top enriched biological pathways.
   - Summarizes the findings into coherent biological hypotheses explaining what is most interesting about the dataset.
```

### Prompt Template 3: Multi-Condition / Screen Analysis (`type = "multicontrast"`)
*Use this prompt when analyzing multiple treatment doses, genotypes, or timepoints (`automating.Rmd`):*

```markdown
You are an expert quantitative proteomics data scientist using `tidyproteomics` (see SKILL.md).

I have a multi-condition experiment with the following pairwise comparisons:
[List comparisons, e.g.:
 - c("knockdown", "control")
 - c("rescue", "control")
 - c("rescue", "knockdown")]

Write a complete R script that:
1. Normalizes the dataset using `select_normalization()`.
2. Automates all two-sample tests in batch using `expression(.pairs = pairs, .method = "limma")`.
3. Automates multi-term GSEA across biological processes and pathways using `enrichment(.pairs = pairs, .terms = c("biological_process", "reactome_pathway"))`.
4. Saves out all figures and tables automatically using `analyze_expressions(labels_column = "gene_name")` and `analyze_enrichments(significance_max = 0.05)`.
5. Generates a comparative expression plot (`plot_compexp()`) to contrast fold changes between conditions.
```

### Prompt Template 4: Venn Diagram & Region Extraction (`type = "extracting"`)
*Use this prompt when investigating shared vs. unique proteomes across conditions (`extracting.Rmd`):*

```markdown
You are an expert quantitative proteomics data scientist using `tidyproteomics` (see SKILL.md).

I have an imported `tidyproteomics` dataset with multiple sample groups:
[List groups, e.g. "control", "drug_low", "drug_high"]

I want to investigate which proteins are uniquely expressed in specific conditions versus shared across groups.

Write a complete R script that:
1. Filters for non-imputed proteins using `subset(imputed == 0)` to ensure true mass spectrometry identifications.
2. Generates and saves a Venn/Euler diagram (`plot_venn(destination = "png")` or `plot_euler(destination = "png")`).
3. Uses `intersection(.include = ..., .exclude = ...)` to extract the underlying data for:
   - Proteins unique to [GROUP_A] (absent in other groups).
   - Proteins shared across all groups.
   - Specific overlap segments of interest.
4. Uses `extract(na.rm = TRUE)` to convert each extracted segment into a tidy long-format tibble and exports them to CSV (`readr::write_csv()`).
5. Prints a concise summary of each extracted segment using `print()`.
```

