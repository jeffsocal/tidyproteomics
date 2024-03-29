% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fun_export_config.R
\name{export_config}
\alias{export_config}
\title{Helper function to export the config file to current project directory}
\usage{
export_config(platform = NULL, analyte = c("proteins", "peptides"))
}
\arguments{
\item{platform}{the source of the data (ProteomeDiscoverer, MaxQuant)}

\item{analyte}{the omics analyte (proteins, peptides)}
}
\value{
success or fail
}
\description{
Helper function to export the config file to current project directory
}
\examples{
library(tidyproteomics)
#\dontrun{
export_config("mzTab", 'peptides')
}

}
