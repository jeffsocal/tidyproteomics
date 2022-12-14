% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/protein_map.R
\name{protein_map}
\alias{protein_map}
\title{Align a peptide data to protein sequences for visualization}
\usage{
protein_map(data = NULL, fasta_path = NULL)
}
\arguments{
\item{data}{a tidyproteomics data-object, specifically of peptide origin}

\item{fasta_path}{a character string representing the path to a fasta file}
}
\value{
a list of protein mappings
}
\description{
Align a peptide data to protein sequences for visualization
}
\examples{
library(dplyr, warn.conflicts = FALSE)
library(tidyproteomics)

hela_protein_map <- hela_peptides \%>\%
   protein_map(fasta = path_to_package_data('fasta'))

}
