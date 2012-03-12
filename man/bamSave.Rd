\name{bamSave}
\alias{bamSave}
\title{bamSave(bamAlign)}
\description{Saves an object of class bamAlignment or bamRange to a BAM-file via an Instance of class bamWriter.}
\usage{bamSave(object, ...)}
\arguments{
  \item{object}{An Instance of class bamWriter}
  \item{\dots}{An Instance of class bamAlign or bamRange}
}
\value{None}
\author{Wolfgang Kaisers}

\examples{
library(rbamtools)
bam_file <- system.file("extdata", "accepted_hits.bam", package="rbamtools")
idx_file <- system.file("extdata", "accepted_hits.bam.bai", package="rbamtools")
reader<-bamReader(bam_file)
isOpen(reader)
align<-getNextAlign(reader)
loadIndex(reader,idx_file)
coords<-as.integer(c(0,0,1874))
index.initialized(reader)
range<-bamRange(reader,coords)
writer<-bamWriter(reader,"test.bam")
bamSave(writer,align)
bamSave(writer,range)
bamClose(reader)
bamClose(writer)
}
\keyword{bamReader}
\keyword{bamAlign}
\keyword{bamRange}
\keyword{bamSave}
