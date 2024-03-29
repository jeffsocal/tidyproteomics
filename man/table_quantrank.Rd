% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_quantrank.R
\name{table_quantrank}
\alias{table_quantrank}
\title{Helper function to quantitation plots}
\usage{
table_quantrank(
  data = NULL,
  accounting = NULL,
  display_filter = c("none", "log2_foldchange", "p_value", "adj_p_value")
)
}
\arguments{
\item{data}{tidyproteomics data object}

\item{accounting}{character string}

\item{display_filter}{a numeric between 0 and 1}
}
\value{
a (tidyproteomics data-object | ggplot-object)
}
\description{
\code{table_quantrank()}
}
\examples{
library(dplyr, warn.conflicts = FALSE)
library(tidyproteomics)
hela_proteins \%>\% plot_quantrank()

hela_proteins \%>\% plot_quantrank(type = 'lines')

hela_proteins \%>\% plot_quantrank(type = 'lines', display_filter = 'log2_foldchange', display_cutoff = 1)

}
