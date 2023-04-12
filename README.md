# tidyproteomics <img src="man/figures/logo.png" align="right" height="139"/>

An R package for the tidy-ing, post processing and analysis of quantitative proteomic data.

Proteomics analysis software, available either through a paid subscription or as an open-source tool, fail to output data in a well conceived tidy format. A majority of these tools generate output formats that have either mixed wide- and long-format data structures, columns headers with messy names and added symbols, and often confusing variable names. This leads researchers to create one-off scripts for cleaning and importing data from various formats, often creating an environment of unmaintained, bespoke code. This package attempts to solve that problem by creating a flexible import tool to unify multiple formats and create an new tidy R object for proteomics analysis.

This package supports at a high level:

-   data importing
-   data filtering
-   data visualization
-   quantitative normalization & imputation
-   two-sample expression & term enrichment analysis
-   protein inference, sequence coverage and visualization

## Data Support

Importing is currently implemented for a few platforms and assume peptide level FDR (at the user's desired level) has already been accounted for. See `vignette("importing")`. Importing is flexible enough to accept other data platforms in flat files (.csv, .tsv, and .xlsx) with a custom configuration.

| Platform           | peptides                        | proteins                  | notes                         |
|-----------------|--------------------|-----------------|-------------------|
| ProteomeDiscoverer | \*.xlsx *peptides export*       | \*.xlsx *proteins export* | requires layout configuration |
| MaxQuant           | evidence.txt                    | proteinGroups.txt         |                               |
| FragPipe           | combined_peptide.tsv            | combined_protein.tsv      |                               |
| Skyline            | \*.csv *MSstats peptide report* |                           | requires MSstats install      |
| DIA-NN             | \*.tsv *peptide report*         |                           |                               |
| mzTab              | \*.mzTab (v1.0.0)               | \*.mzTab (v1.0.0)         | does not track MBR            |

## Ease of Use

This package supports the same [syntactic sugar](https://en.wikipedia.org/wiki/Syntactic_sugar) utilized in the tidy-verse functions like filter, and introduces the `%like%` operator, see `vignette("subsetting")` . These operations can extend to all aspects of the data set, including sample names, protein IDs, annotations and accountings like *match_between_runs* and *num_peptides*.

| operator | description          | example                                          |
|-----------------|-----------------|--------------------------------------|
| ==       | equals               | `sample == 'wt'` , `match_between_runs == FALSE` |
| !=       | does not equal       | `biological_function != 'DNA metabolism'`        |
| \<, \>   | less, greater than   | `num_unique_peptides >= 2`                       |
| %like%   | contains             | `description %like% 'ribosome'`                  |
| ! %like% | does not contain     | `!description %like% 'ribosome'`                 |
| ---      | --- *expression* --- | ---                                              |
| /        | ratio                | `experiment / control`                           |

Expression analysis also utilizes this type of syntax when referencing samples for analysis. For example `data %>% expression(knockdown/control)` would know to run the differential expression of the sample *ko* with respect to the sample *wt* such that positive log2 difference would be up-expressed in *ko* and a negative log2 differences would be down-expressed in *ko*.

## Installation

To install, open R and type:

``` r
install.packages("devtools")
devtools::install_github("jeffsocal/tidyproteomics")
```

You will also need the Bioconductor packages [limma](https://bioconductor.org/packages/release/bioc/html/limma.html), [qvalue](https://bioconductor.org/packages/release/bioc/html/qvalue.html), and [fgsea](https://bioconductor.org/packages/release/bioc/html/fgsea.html), to install these type:

``` r
install.packages("BiocManager")
BiocManager::install(c("limma","qvalue","fgsea"))
```

**NOTE:** There are several other packages required that will be prompted and automatically downloaded from CRAN when installing. Depending on your current system some packages require the installation of OS level libraries for advanced math and string manipulation.

## Get Started
If you want to test the features and capabilities of this project and get a better understanding of how it works, you can experiment with the two datasets that are included within the package. These datasets have been specially selected to provide you with the necessary working examples that will help you to get started with the project.

The first dataset is a protein-level dataset, and the second one is a small peptide-level dataset. Both datasets are easily accessible after loading the package. You can use these datasets to practice and explore the different features and functionalities of the package. They are also used throughout the vignettes and code documentation to provide you with clear and concise examples.

Protein level data.
```r
library(tidyproteomics)

hela_protiens

── Quantitative Proteomics Data Object ──

Origin          ProteomeDiscoverer 
                proteins (10.76 MB) 
Composition     6 files 
                2 samples (control, knockdown) 
Quantitation    7055 proteins 
                4 log10 dynamic range 
                28.8% missing values 
Accounting      (4) match_between_runs num_peptides num_unique_peptides num_psms 
Annotations     (9) description gene_id_entrez gene_id_ensemble gene_name biological_process
                cellular_component molecular_function wiki_pathway
                reactome_pathway 
                
```

Peptide level data.
```r
library(tidyproteomics)

hela_peptides
── Quantitative Proteomics Data Object ──

Origin          ProteomeDiscoverer 
                peptides (188.87 kB) 
Composition     6 files 
                2 samples (control, knockdown) 
Quantitation    8 proteins 
                3.3 log10 dynamic range 
                21.3% missing values 
Accounting      (3) modifications match_between_runs num_psms 
```


## Import Your Data
Its simple to get started. Make a new project, drop your raw data in a folder labeled *data*. For more information see `vignette("workflow-simple")`

``` r
library(tidyproteomics)

data <- "./data/some_ProteomeDiscoverer_data.xlsx" %>%
  import("ProteomeDiscoverer", "proteins")
  
data %>% summary("samples")
```


