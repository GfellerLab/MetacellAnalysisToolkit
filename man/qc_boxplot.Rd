% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/qc_boxplot.R
\name{qc_boxplot}
\alias{qc_boxplot}
\title{Generate boxplot to visualize QC metrics distributions.}
\usage{
qc_boxplot(mc.obj, qc.metrics = NULL, split.by = "orig.ident", y.lim = NULL)
}
\arguments{
\item{mc.obj}{A Seurat object containing the metacells data.}

\item{qc.metrics}{Vector of strings indicating which QC metric should be considered in the mc.obj metadata.}

\item{split.by}{(optional): String indicating if the boxplot should be splitted based on a metacell annotation available in the mc.obj metadata dataframe.}

\item{y.lim}{(optional): Y axis limit. 
By default, the orig.ident variable is used.}
}
\value{
Boxplots representing the distribution of the qc.metrics.
}
\description{
\code{qc_boxplot} 
This function generate boxplots to visualize QC metrics distributions across metacells.
}
\examples{
qc_boxplot(mc.obj = CD34_mc, qc.metrics = c("size", "celltype_purity"))
}
