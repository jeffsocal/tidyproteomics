% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/analysis_enrichment_gsea.R
\name{enrichment_gsea}
\alias{enrichment_gsea}
\title{A function for evaluating term enrichment via GSEA}
\usage{
enrichment_gsea(data_expression = NULL, data = NULL, term_group = NULL)
}
\arguments{
\item{data_expression}{a tibble from and two sample expression difference analysis}

\item{data}{tidyproteomics data object}

\item{term_group}{a character string referencing "term" in the annotations table}
}
\value{
a tibble
}
\description{
A function for evaluating term enrichment via GSEA
}
