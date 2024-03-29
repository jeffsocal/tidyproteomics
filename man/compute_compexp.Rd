% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_compexp.R
\name{compute_compexp}
\alias{compute_compexp}
\title{Helper function to analysis between two expression tests}
\usage{
compute_compexp(
  table_a = NULL,
  table_b = NULL,
  log2fc_min = 2,
  log2fc_column = "log2_foldchange",
  significance_max = 0.05,
  significance_column = "adj_p_value",
  labels_column = "protein"
)
}
\arguments{
\item{table_a}{a tibble}

\item{table_b}{a tibble}

\item{log2fc_min}{a numeric defining the minimum log2 foldchange to highlight.}

\item{log2fc_column}{a character defining the column name of the log2 foldchange values.}

\item{significance_max}{a numeric defining the maximum statistical significance to highlight.}

\item{significance_column}{a character defining the column name of the statistical significance values.}

\item{labels_column}{a character defining the column name of the column for labeling.}
}
\value{
a list
}
\description{
Helper function to analysis between two expression tests
}
