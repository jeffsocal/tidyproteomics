## CRAN Submission Comments

### Package: tidyproteomics 1.9.15

### Test environments
* local macOS (R 4.x)
* win-builder (R-devel, x86_64-w64-mingw32, Windows Server 2022)

### R CMD check results
Status: 0 ERRORs, 0 WARNINGs, 2 NOTEs

### Explanation of NOTEs

* **checking CRAN incoming feasibility ... NOTE**
  - **New submission**: This is the initial submission of `tidyproteomics` to CRAN.
  - **Possibly misspelled words in DESCRIPTION**:
    `Proteomic`, `proteomic`, `spectrometry`.
    These are correct scientific and analytical terminology (mass spectrometry, quantitative proteomics) and are listed in `inst/WORDLIST`.
  - **Availability using Additional_repositories specification**:
    Bioconductor repository (`https://bioconductor.org/packages/release/bioc`) is specified in `DESCRIPTION` for Bioconductor dependencies (`limma`, `fgsea`).
  - **Dependencies**:
    Optional packages (`BiocParallel`, `e1071`, `eulerr`, `pheatmap`, `randomForest`, `RColorBrewer`, `readxl`, `writexl`) are placed in `Suggests` and handled conditionally via `rlang::check_installed()`, keeping non-default `Imports` at 20.

* **checking examples ... NOTE**
  - **Long-running examples**:
    Example runtimes have been addressed:
    - Long-running examples in `plot_protein()`, `protein_map()`, `plot_counts()`, and `impute()` have been wrapped in `\donttest{}`.
    - The example in `summary.tidyproteomics()` was updated to use a subset (`hela_proteins %>% subset(description %like% 'ribosome')`) to ensure minimal compute time during checks.
