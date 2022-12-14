% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data_annotate.R
\name{annotate}
\alias{annotate}
\title{Main function for adding annotations to a tidyproteomics data-object}
\usage{
annotate(
  data = NULL,
  annotations = NULL,
  duplicates = c("replace", "merge", "leave")
)
}
\arguments{
\item{data}{a tidyproteomics data list-object}

\item{annotations}{a character string vector}

\item{duplicates}{a character string, how to handle duplicate terms}
}
\value{
a tidyproteomics data list-object
}
\description{
Main function for adding annotations to a tidyproteomics data-object
}
