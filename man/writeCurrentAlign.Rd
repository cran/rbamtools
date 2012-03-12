\name{writeCurrentAlign}
\alias{writeCurrentAlign}
\alias{insertPreCurrent}
\alias{insertPastCurrent}
\title{writeCurrentAlign, insertPreCurrent, insert PastCurrent: Writing bamAlign data to bamRange list.}
\description{The three functions can be used to write bamAlign data to bamRange object.
The function writeCurrentAlign overwrites the data of the current bamRange align. When no current align is set (i.e. ==NULL), nothing is written.

The functions insertPreCurrent and insertPastCurrent insert an bamAlign object into the bamRange list. When no current align is set, insertPreCurrent inserts the object as the last object and insertPastCurrent inserts the object as the first align object in the bamRange list.
}
\usage{
writeCurrentAlign(object, value)
insertPreCurrent(object,value)
insertPastCurrent(object,value)
}
\arguments{
  \item{object}{An instance of class bamRange.}
  \item{value}{An instance of class bamAlign.}
}
\value{None.}
\author{Wolfgang Kaisers}
\seealso{bamRange}
\examples{
bam_file <- system.file("extdata", "accepted_hits.bam", package="rbamtools")
idx_file <- system.file("extdata", "accepted_hits.bam.bai", package="rbamtools")
reader<-bamReader(bam_file)
isOpen(reader)
align<-getNextAlign(reader)
loadIndex(reader,idx_file)
index.initialized(reader)
coords<-as.integer(c(0,0,1874))
range<-bamRange(reader,coords)
align<-getNextAlign(range)
insertPastCurrent(range,align)
insertPreCurrent(range,align)
writeCurrentAlign(range,align)
bamClose(reader)
}
\keyword{bamRange}
\keyword{bamAlign}
