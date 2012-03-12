\name{bamRange}
\alias{bamRange}
\title{bamRange(bamReader, coordinates)}
\description{bamRange retrieves all alignments in the opened}
\usage{bamRange(reader, coords)}
\arguments{
 \item{reader}{An instance of bamReader. Must be opened and contain initialized index}
  \item{coords}{Integer vector of length 3: coords=c(refid,start,stop)}
}
\details{The method returns a list of bamAlign's from which overlap with the specified region.}
\value{An instance of class bamRange which can be accessed sequentially, modified or written to a BAM-file.}
\examples{
library(rbamtools)
bam_file <- system.file("extdata", "accepted_hits.bam", package="rbamtools")
idx_file <- system.file("extdata", "accepted_hits.bam.bai", package="rbamtools")
reader<-bamReader(bam_file)
isOpen(reader)
align<-getNextAlign(reader)
loadIndex(reader,idx_file)
index.initialized(reader)
#################################################
# Find the appropriate refID
rdf<-getRefData(reader)
head(rdf)
#################################################
# The sequence length (LN) also determines valid
# range for start and stop coordinates
# Invalid refid-, start- or stop-coordinates will
# release an error.
# coords: refid=0, start=0, stop=1874
coords<-as.integer(c(0,0,1874))
range<-bamRange(reader,coords)
align<-getNextAlign(range)
bamClose(reader)
}
\keyword{bamRange}
\keyword{index}
