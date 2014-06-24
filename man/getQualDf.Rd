\name{getQualDf}
\alias{getQualDf}
\alias{getQualDf-methods}
\alias{getQualDf,bamRange-method}
\alias{getQualQuantiles}
\alias{getQualQuantiles-methods}
\alias{getQualQuantiles,bamRange-method}
\alias{plotQualQuant}
\alias{plotQualQuant-methods}
\alias{plotQualQuant,bamRange-method}
\title{Read and display Phred qualities from bamRange}
\description{\code{getQualDf} takes a \code{bamRange} and returns a data.frame (128 rows, number of columns=length of the longest sequence in range).\code{getQualDf} counts occurrences for every sequence position (column) and every phred value (row). \code{getQualQuantiles} takes a \code{bamReader} and a vector of quantiles (must be between 0 and 1) and returns a data.frame. The data.frame contains one row for each quantile and also as many columns as the maximum sequence length. \code{plotQualQuant} plots the values for quanties 0.1,0.25,0.5,0.75 and 0.9.}
\usage{getQualDf(object,prob=FALSE,...)}
\arguments{
  \item{object}{bamRange.}
  \item{prob}{logical. When TRUE each column is divided by its sum. The column sums are added as names attribute \code{col.sum}}
  \item{...}{(currently unused)}
}
\value{data.frame}
\details{Phred values are truncated by 127 (the maximum which can be represented by ASCII values). The function runs down each colum (sequence position) and returns the row index where the quantile exceeds the cumlated column values.}
\author{Wolfgang Kaisers}
\examples{
# A) Read bamRange
bam<-system.file("extdata","accepted_hits.bam",package="rbamtools")
reader<-bamReader(bam,idx=TRUE)
coords<-as.integer(c(0,0,249250621))
range<-bamRange(reader,coords)
bamClose(reader)
# B) getQualDf
qdf<-getQualDf(range)
qdf[32:38,1:15]
qdr<-getQualDf(range,prob=TRUE)
# C) getQualQuantiles
quantiles<-c(0.1,0.25,0.5,0.75,0.9)
qt<-getQualQuantiles(range,quantiles)
# D) Plot
plotQualQuant(range)
}
\keyword{getQualDf}
\keyword{bamRange}