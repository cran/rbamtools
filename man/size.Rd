\name{size}
\alias{size}
\alias{size-methods}
\title{Returns size of bamRange list.}
\description{Returns number of bamAligns stored in bamRange.}
\usage{size(object)}
\arguments{
  \item{object}{bamRange object}
}
\value{Integer.}
\references{
\verb{The SAM Format Specification (v1.4-r985). The SAM Format Specification Working Group. September 7, 2011.}
\url{http://www.samtools.sourceforge.net/SAM1.pdf}
}
\author{Wolfgang Kaisers}
\examples{
library(rbamtools)
bam_file <- system.file("extdata", "accepted_hits.bam", package="rbamtools")
idx_file <- system.file("extdata", "accepted_hits.bam.bai", package="rbamtools")
reader<-bamReader(bam_file)
isOpen(reader)
align<-getNextAlign(reader)
loadIndex(reader,idx_file)
index.initialized(reader)
coords<-as.integer(c(0,0,1874))
range<-bamRange(reader,coords)
size(range)
bamClose(reader)
}
\keyword{bamRange}
\keyword{size}
