% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fun_data_io.R
\name{read_data}
\alias{read_data}
\title{Read data by format type}
\usage{
read_data(path = NULL, ...)
}
\arguments{
\item{path}{a path character string}

\item{...}{readr passthrough options}
}
\value{
tibble
}
\description{
\code{read_data()} is a helper function that assumes the format type of the data
table by checking the ending of path string
}
