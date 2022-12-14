% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fun_hdf.R
\name{hdf}
\alias{hdf}
\title{Helper function to take the head of a tibble and display as a data.frame}
\usage{
hdf(x, n = 5)
}
\arguments{
\item{x}{a tibble}

\item{n}{display up to the nth row}
}
\value{
a data frame
}
\description{
Helper function to take the head of a tibble and display as a data.frame
}
\examples{
library(tidyproteomics)
x <- tibble::tibble(a = 1:10, b = 11:20)
hdf(x)
hdf(x, n = 3)
}
