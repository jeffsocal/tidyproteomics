# tidyproteomics <a href=''><img src='man/figures/logo.png' align="right" height="139" /></a>

An R package for the post processing and analysis of quantitative
proteomic data - effectively wiring inputs and outputs for complex
analyses akin to bundled copper wires that are used to create a circuit
for residential electrical services. This package supports at a high
level:

-   data filtering
-   data visualization
-   quantitative normalization & imputation
-   two-sample expression & term enrichment analysis

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
`vignette("workflow-simple")`

``` r
library(tidyproteomics)

data <- "./data/some_ProteomeDiscoverer_data.xlsx" %>%
  import("ProteomeDiscoverer", "proteins")
  
data %>% summary("samples")
```
