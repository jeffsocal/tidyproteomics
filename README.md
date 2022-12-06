# tidyproteomics <a href=''><img src="man/figures/logo.png" align="right" height="139"/></a>

An R package for the post processing and analysis of quantitative
proteomic data. In reality, this is just a package of functions and
plots that I use every day for the analysis of quantitative proteomic
data. Currently inputs are implemented for ProteomeDiscoverer, MaxQuant,
Skyline and DIA-NN, and assume peptide level FDR has already been
accounted for. This package supports at a high level:

-   data filtering
-   data visualization
-   quantitative normalization & imputation
-   two-sample expression & term enrichment analysis
-   protein inference, sequence coverage and visualization

This package supports the same [syntactic
sugar](https://en.wikipedia.org/wiki/Syntactic_sugar) utilized in the
tidy-verse functions like filter, and introduces the `%like%` operator,
see `vignette("subsetting")` . These operations can extend to all
aspects of the data set, including sample names, protein IDs,
annotations and accountings like *match_between_runs* and
*num_peptides*.

| operator        | description          | example                                          |
|-------------------|-------------------|-----------------------------------|
| ==              | equals               | `sample == 'wt'` , `match_between_runs == FALSE` |
| !=              | does not equal       | `biological_function != 'DNA metabolism'`        |
| \<, \>          | less, greater than   | `num_unique_peptides >= 2`                       |
| %like%, %!like% | contains             | `description %like% 'ribosome'`                  |
| ---             | --- *expression* --- | ---                                              |
| /               | ratio                | `experiment / control`                           |

Expression analysis also utilizes this type of syntax when referencing
samples for analysis. For example `data %>% expression(knockdown/control)` would
know to run the differential expression of the sample *ko* with respect
to the sample *wt* such that positive log2 difference would be
up-expressed in *ko* and a negative log2 differences would be
down-expressed in *ko*.

## Installation

To install, open R and type:

``` r
install.packages("devtools")
devtools::install_github("jeffsocal/tidyproteomics")
devtools::install_github("jeffsocal/rfasta")
```

You will also need the Bioconductor packages
[limma](https://bioconductor.org/packages/release/bioc/html/limma.html),
[qvalue](https://bioconductor.org/packages/release/bioc/html/qvalue.html),
and
[fgsea](https://bioconductor.org/packages/release/bioc/html/fgsea.html),
to install these type:

``` r
install.packages("BiocManager")
BiocManager::install(c("limma","qvalue","fgsea"))
```

## Get Started

Its simple to get started. Make a new project, drop your raw data in a
folder labeled *data*. For more information see
`vignette("01-workflow-simple")`

``` r
library(tidyproteomics)

data <- "./data/some_ProteomeDiscoverer_data.xlsx" %>%
  import("ProteomeDiscoverer", "proteins")
  
data %>% summary("samples")
```
