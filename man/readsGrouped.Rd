% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CNV_Plot.R
\name{readsGrouped}
\alias{readsGrouped}
\title{Title: Grouping spots based on the number of reads.}
\usage{
readsGrouped(data1, data2, Treads = 8000, NumSpt = 8, locpath)
}
\arguments{
\item{data1}{Original filtered_feature_bc.csv}

\item{data2}{grouped_spotSummary.csv}

\item{Treads}{Minimum number of reads in a new spot}

\item{NumSpt}{Number of spots used for grouping}

\item{locpath}{The tissue_positions_list.csv}
}
\value{
Grouped spots files as new input
}
\description{
Title: Grouping spots based on the number of reads.
}
