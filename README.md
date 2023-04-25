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

+--------------------+---------------------------------+---------------------------+-------------------------------+
| Platform           | peptides                        | proteins                  | notes                         |
+====================+=================================+===========================+===============================+
| ProteomeDiscoverer | \*.xlsx *peptides export*       | \*.xlsx *proteins export* | requires layout configuration |
+--------------------+---------------------------------+---------------------------+-------------------------------+
| MaxQuant           | evidence.txt                    | proteinGroups.txt         |                               |
+--------------------+---------------------------------+---------------------------+-------------------------------+
| FragPipe           | combined_peptide.tsv            | combined_protein.tsv      |                               |
+--------------------+---------------------------------+---------------------------+-------------------------------+
| Skyline            | \*.csv *MSstats peptide report* |                           | requires MSstats install      |
+--------------------+---------------------------------+---------------------------+-------------------------------+
| DIA-NN             | \*.tsv *peptide report*         |                           |                               |
+--------------------+---------------------------------+---------------------------+-------------------------------+
| mzTab              | \*.mzTab (v1.0.0)               | \*.mzTab (v1.0.0)         | does not track MBR            |
+--------------------+---------------------------------+---------------------------+-------------------------------+

## Ease of Use

This package supports the same [syntactic sugar](https://en.wikipedia.org/wiki/Syntactic_sugar) utilized in the tidy-verse functions like filter, and introduces the `%like%` operator, see `vignette("subsetting")` . These operations can extend to all aspects of the data set, including sample names, protein IDs, annotations and accountings like *match_between_runs* and *num_peptides*.

+--------------------------------------------------------------------------+----------------------+--------------------------------------------------------------------------+
| operator                                                                 | description          | example                                                                  |
+==========================================================================+======================+==========================================================================+
| ==                                                                       | equals               | `sample == 'wt'` , `match_between_runs == FALSE`                         |
+--------------------------------------------------------------------------+----------------------+--------------------------------------------------------------------------+
| !=                                                                       | does not equal       | `biological_function != 'DNA metabolism'`                                |
+--------------------------------------------------------------------------+----------------------+--------------------------------------------------------------------------+
| \<, \>                                                                   | less, greater than   | `num_unique_peptides >= 2`                                               |
+--------------------------------------------------------------------------+----------------------+--------------------------------------------------------------------------+
| %like%                                                                   | contains             | `description %like% 'ribosome'`                                          |
+--------------------------------------------------------------------------+----------------------+--------------------------------------------------------------------------+
| ! %like%                                                                 | does not contain     | `!description %like% 'ribosome'`                                         |
+--------------------------------------------------------------------------+----------------------+--------------------------------------------------------------------------+
| /                                                                        | ratio _(expression)_ | `experiment / control`                                                   |
+--------------------------------------------------------------------------+----------------------+--------------------------------------------------------------------------+

Expression analysis also utilizes this type of syntax when referencing samples for analysis. For example `data %>% expression(knockdown/control)` would know to run the differential expression of the sample *ko* with respect to the sample *wt* such that positive log2 difference would be up-expressed in *ko* and a negative log2 differences would be down-expressed in *ko*.

## Nomenclature

The tidyprotoemics package uses a conserved naming convention that facilitates a "contract" with funtions that operate on the data, making it easier to maintain and modify the code. This naming convention ensures that functions are generic to their inputs, being able to operate on a variety of data without having to explicitly name each column and variable which may be named different between data sets. All of the effort of standardizing the variable names in the tidyproteomics data object is done upfront when importing the data. See `vignette("importing")`.

This package divides data into tables similar to those in SQL, making it easy to organize and analyze data. This table structure offers a convenient way to explore relationships between different variables and to filter data according to specific criteria.

### Experimental Table

+-------------+-----------------------------------------------------------------------------------------------------------------------------------------+
| Variable    | Reference                                                                                                                               |
+=============+=========================================================================================================================================+
| sample_id   | A checksum id on the individual data imported, allows for an unambiguous identifier usefull in differentiating samples of the same name |
+-------------+-----------------------------------------------------------------------------------------------------------------------------------------+
| import_file | The file referenced at import, allows for multiple files to be merged                                                                   |
+-------------+-----------------------------------------------------------------------------------------------------------------------------------------+
| sample_file | The individual MS file, or for TMT the channel, that contains data for a single observation of a given protein                          |
+-------------+-----------------------------------------------------------------------------------------------------------------------------------------+
| sample      | The name of the MS data collected for a given class, collection, of a specimen or a cell culture, etc.                                  |
+-------------+-----------------------------------------------------------------------------------------------------------------------------------------+
| replicate   | A tidyproteomics derived integer assignment to a collection of same-named samples                                                       |
+-------------+-----------------------------------------------------------------------------------------------------------------------------------------+

[Example]{.underline}

``` r
library(tidyproteomics)

hela_proteins$experiments

# A tibble: 6 × 5
  sample_id import_file                sample_file sample    replicate
  <chr>     <chr>                      <chr>       <chr>     <chr>    
1 e9b20ea7  p97KD_HCT116_proteins.xlsx f1          control   1        
2 ef79cc4c  p97KD_HCT116_proteins.xlsx f4          control   2        
3 eebba67b  p97KD_HCT116_proteins.xlsx f5          control   3        
4 ebf4b0fe  p97KD_HCT116_proteins.xlsx f2          knockdown 1        
5 ea36dac9  p97KD_HCT116_proteins.xlsx f3          knockdown 2        
6 ecfd1822  p97KD_HCT116_proteins.xlsx f6          knockdown 3   
```

### Quantitation Table

+--------------------------------+-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| Variable                       | Reference                                                                                                                                                                                                                                                                                                                   |
+================================+=============================================================================================================================================================================================================================================================================================================================+
| sample_id                      | A checksum id on the individual data imported, allows for an unambiguous identifier useful in differentiating samples of the same name                                                                                                                                                                                      |
+--------------------------------+-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| sample                         | The name of the MS data collected for a given class, collection of sampling of tissue, a culture etc.                                                                                                                                                                                                                       |
+--------------------------------+-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| replicate                      | A tidyproteomics derived integer assignment to a collection of same-named samples                                                                                                                                                                                                                                           |
+--------------------------------+-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| *identifier*                   | The identifier of the "thing" being measured.                                                                                                                                                                                                                                                                               |
|                                |                                                                                                                                                                                                                                                                                                                             |
| protein                        | In the case of protein-level data, this is simply `protein` and is usually populated with UniProt accession numbers, names, or other protein identifiers.                                                                                                                                                                   |
|                                |                                                                                                                                                                                                                                                                                                                             |
| protein, peptide, modification | In the case of peptide-level data, this is a multiple of `protein` (accession number), `peptide` (amino acid sequence) and `modification` (post-translational modification).                                                                                                                                                |
+--------------------------------+-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| abundance_raw                  | The imported quantitative measure for each protein or peptide. The tidyproteomics data object stores abundance values for the raw data and each normalization. This allows for the direct comparison of normalization methods against the raw data, and for one set of abundance values to be used for subsequent analyses. |
+--------------------------------+-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+

[Example]{.underline}

``` r
library(tidyproteomics)

hela_proteins$quantitation

# A tibble: 42,330 × 5
   sample_id sample    replicate protein abundance_raw
   <chr>     <chr>     <chr>     <chr>           <dbl>
 1 e9b20ea7  control   1         Q15149    1011259992.
 2 ef79cc4c  control   2         Q15149    1093277593.
 3 eebba67b  control   3         Q15149     980809516.
 4 ebf4b0fe  knockdown 1         Q15149    1410445367.
 5 ea36dac9  knockdown 2         Q15149    1072305561.
 6 ecfd1822  knockdown 3         Q15149    1486561518.
# … with 42,324 more rows
```

### Annotation Table

+--------------+------------------------------------------------------------------------------+
| Variable     | Reference                                                                    |
+==============+==============================================================================+
| *identifier* | The identifier of the "thing" being measured.                                |
+--------------+------------------------------------------------------------------------------+
| term         | The name of the collection of annotations, such as gene_name or description. |
+--------------+------------------------------------------------------------------------------+
| annotation   | The descriptive annotation.                                                  |
+--------------+------------------------------------------------------------------------------+

[Example]{.underline}

``` r
library(tidyproteomics)

print(hela_proteins$annotations)

# A tibble: 63,495 × 3
   protein term               annotation                                                                                                          
   <chr>   <chr>              <chr>                                                                                                               
 1 Q15149  description        Plectin OS=Homo sapiens OX=9606 GN=PLEC PE=1 SV=3                                                                   
 2 Q15149  gene_id_entrez     5339                                                                                                                
 3 Q15149  gene_id_ensemble   ENSG00000178209                                                                                                     
 4 Q15149  gene_name          PLEC                                                                                                                
 5 Q15149  biological_process cell differentiation;cell growth;cellular homeostasis;coagulation;defense response;metabolic process
# … with 63,490 more rows
```

### 

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

``` r
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

``` r
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
