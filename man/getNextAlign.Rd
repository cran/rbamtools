\name{getNextAlign}
\alias{getNextAlign}
\title{
Reads next Align from bamReader or bamRange.
}
\description{
An open bamReader and a bamRange maintain a pointer to a current align. getNextAlign shifts the current align to the next following align and returns the value of the current align.

When there is no align past the current align (e.g. EOF), the getNextAling returns NULL. When a corrupted BAM-file is detected, an error is thrown.
}
\usage{getNextAlign(object)}
\arguments{
  \item{object}{An object of class bamReader or bamRange.}
}
\value{An object of class bamAlign or NULL.}
\author{Wolfgang Kaisers}
\examples{
library(rbamtools)
bam_file <- system.file("extdata", "accepted_hits.bam", package="rbamtools")
idx_file <- system.file("extdata", "accepted_hits.bam.bai", package="rbamtools")
reader<-bamReader(bam_file)
isOpen(reader)
align<-getNextAlign(reader)
if(is.null(align))
	print("End of File!\n")
loadIndex(reader,idx_file)
index.initialized(reader)
coords<-as.integer(c(0,0,1874))
range<-bamRange(reader,coords)
align<-getNextAlign(range)

align<-getNextAlign(range)
bamClose(reader)
}
\keyword{bamAlign}
\keyword{bamReader}
\keyword{bamRange}
