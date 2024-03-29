% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_proportion.R
\name{plot_proportion}
\alias{plot_proportion}
\title{Plot proportional expression values}
\usage{
plot_proportion(
  data = NULL,
  ...,
  log2fc_column = "log2_foldchange",
  log2fc_min = 2,
  significance_column = "adj_p_value",
  significance_max = 0.05,
  proportion_column = "proportional_expression",
  proportion_min = 0.01,
  labels_column = NULL,
  label_significance = TRUE,
  show_pannels = FALSE,
  show_lines = TRUE,
  show_fc_scale = TRUE,
  point_size = NULL,
  color_positive = "dodgerblue",
  color_negative = "firebrick1",
  destination = "plot",
  height = 5,
  width = 8
)
}
\arguments{
\item{data}{a tidyproteomics data object}

\item{...}{two sample comparison}

\item{log2fc_column}{a character defining the column name of the log2 foldchange values.}

\item{log2fc_min}{a numeric defining the minimum log2 foldchange to highlight.}

\item{significance_column}{a character defining the column name of the statistical significance values.}

\item{significance_max}{a numeric defining the maximum statistical significance to highlight.}

\item{proportion_column}{a character defining the column name of the proportional expression values.}

\item{proportion_min}{a numeric defining the minimum proportional expression to highlight.}

\item{labels_column}{a character defining the column name of the column for labeling.}

\item{label_significance}{a boolean for labeling values below the significance threshold.}

\item{show_pannels}{a boolean for showing colored up/down expression panels.}

\item{show_lines}{a boolean for showing threshold lines.}

\item{show_fc_scale}{a boolean for showing the secondary foldchange scale.}

\item{point_size}{a numeric for shanging the point size.}

\item{color_positive}{a character defining the color for positive (up) expression.}

\item{color_negative}{a character defining the color for negative (down) expression.}

\item{destination}{a character string}

\item{height}{a numeric}

\item{width}{a numeric}
}
\value{
a ggplot2 object
}
\description{
\code{plot_proportion()} is a GGplot2 implementation for plotting the expression differences
as foldchange ~ scaled abundance. This allows for the visualization of selected
proteins See also \code{plot_volcano()}. This function can
take either a tidyproteomics data object or a table with the required headers.
}
\examples{
library(dplyr, warn.conflicts = FALSE)
library(tidyproteomics)
hela_proteins \%>\%
   expression(knockdown/control) \%>\%
   plot_proportion(knockdown/control, log2fc_min = 0.5, significance_column = 'p_value')

# generates the same out come
# hela_proteins \%>\%
#    expression(knockdown/control) \%>\%
#    export_analysis(knockdown/control, .analysis = 'expression) \%>\%
#    plot_proportion(log2fc_min = 0.5, significance_column = 'p_value')

# display the gene name instead
hela_proteins \%>\%
   expression(knockdown/control) \%>\%
   plot_proportion(knockdown/control, log2fc_min = 0.5, significance_column = 'p_value', labels_column = "gene_name")

}
