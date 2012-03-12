\name{getRefCount}
\alias{getRefCount}
\title{getRefCount: Returns number of reference sequences.}
\description{Returns number of reference sequences as given in BAM-header.}
\usage{getRefCount(object)}
\arguments{
  \item{object}{An object of class bamReader}
}
\value{Integer}
\references{
\verb{The SAM Format Specification (v1.4-r985). The SAM Format Specification Working Group. September 7, 2011.}
\url{http://www.samtools.sourceforge.net/SAM1.pdf}
}
\author{Wolfgang Kaisers}
\examples{
library(rbamtools)
bam_file <- system.file("extdata", "accepted_hits.bam", package="rbamtools")
reader<-bamReader(bam_file)
if(isOpen(reader))
	getRefCount(reader)
bamClose(reader)
}
\keyword{bamReader}
