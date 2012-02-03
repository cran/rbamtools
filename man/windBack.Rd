\name{windBack}
\alias{windBack}
\title{windBack sets the current align of bamRange before the first align (=past the last align) in the bamRange align list.}
\description{
A bamRange object maintains a double-linked list of aligns and a pointer to a current align. windBack shifts the current align before the first align (=past the last align). When there is no align, nothing happens.
}
\usage{
windBack(object)
}
\arguments{
  \item{object}{An object of class bamRange.}
}
\value{NULL.}
\author{Wolfgang Kaisers}
\examples{
library(rbamtools)
bam_file <- system.file("extdata", "accepted_hits.bam", package="rbamtools")
idx_file <- system.file("extdata", "accepted_hits.bam.bai", package="rbamtools")
reader<-bamReader(bam_file)
isOpen(reader)
loadIndex(reader,idx_file)
index.initialized(reader)
coords<-as.integer(c(0,0,1874))
range<-bamRange(reader,coords)
align<-getNextAlign(range)
align<-windBack(range)
align<-getNextAlign(range)
bamClose(reader)
}
\keyword{bamAlign}
\keyword{bamReader}
\keyword{bamRange}
