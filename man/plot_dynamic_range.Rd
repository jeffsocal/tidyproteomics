% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_dynamic_range.R
\name{plot_dynamic_range}
\alias{plot_dynamic_range}
\title{Plot CVs by abundance}
\usage{
plot_dynamic_range(data = NULL, ...)
}
\arguments{
\item{data}{tidyproteomics data object}

\item{...}{passthrough for ggsave see \code{plotting}}
}
\value{
a (tidyproteomics data-object | ggplot-object)
}
\description{
\code{plot_dynamic_range()} is a GGplot2 implementation for plotting the normalization
effects on CVs by abundance, visualized as a 2d density plot. Layered on top
is a loess smoothed regression of the CVs by abundance, with the median CV
shown in \emph{red} and the dynamic range represented as a box plot on top. The
point of this plot is to examine how CVs were minimized through out the abundance
profile. Some normalization methods function well at high abundance yet leave
retain high CVs at lower abundance.
}
\examples{
library(dplyr, warn.conflicts = FALSE)
library(tidyproteomics)
hela_proteins \%>\%
  normalize(.method = c("linear", "loess", "randomforest")) \%>\%
  plot_dynamic_range()

}
