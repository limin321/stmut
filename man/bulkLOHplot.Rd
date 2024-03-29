% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/allelicImbalance.R
\name{bulkLOHplot}
\alias{bulkLOHplot}
\title{Allelic Imbalance or loss of heterozygocity plot.}
\usage{
bulkLOHplot(centmere, alle_imbal)
}
\arguments{
\item{centmere}{hg38 centromere data}

\item{alle_imbal}{tumorshif/shift, which is the variant allele frequency deviated from 50\% VAF}
}
\value{
An arm-level LOH plot
}
\description{
Allelic Imbalance or loss of heterozygocity plot.
}
\examples{
\dontrun{
centm <- read.csv(system.file("extdata/", "hg38_centromereSimple.bed", package = "stmut"),
  sep = "\t", header = FALSE
)
loh <- read.table(system.file("extdata/", "p6LOH.txt", package = "stmut"),
  sep = "\t",
  header = TRUE
)
loh <- loh[, c(1, 2, 3, 10)]
colnames(loh) <- c("chromosome", "start", "end", "tumorshift")
loh$chromosome <- substr(loh$chromosome, 4, nchar(loh$chromosome))
bulkLOHplot(centmere = centm, alle_imbal = loh)
}
}
