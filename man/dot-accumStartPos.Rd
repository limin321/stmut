% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/allelicImbalance.R
\name{.accumStartPos}
\alias{.accumStartPos}
\title{Set genomic start position to accumulated start position.}
\usage{
.accumStartPos(dataf)
}
\arguments{
\item{dataf}{the df must have 2 cols, "chromosome" character show as 1,2,...,"X","Y",
the "start" column must be numeric.}
}
\value{
A data with accumulated start position, called in bulkLOHplot func.
}
\description{
Set genomic start position to accumulated start position.
}
\examples{
\dontrun{
loh <- read.table(system.file("extdata/", "p6LOH.txt", package = "stmut"),
  sep = "\t",
  header = TRUE
)
loh <- loh[, c(1, 2, 3, 10)]
colnames(loh) <- c("chromosome", "start", "end", "tumorshift")
loh$chromosome <- substr(loh$chromosome, 4, nchar(loh$chromosome))
df8 <- accumStartPos(loh)
}
}
