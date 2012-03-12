\name{getPrevAlign}
\alias{getPrevAlign}
\title{Reads previous Align from bamRange.}
\description{
A bamRange object maintains a pointer to a current align. getPrevAlign shifts the current align to the previous align and then returns the value of the current align. When there is no align before the current align (e.g. first align in list), the getPrevAlign returns NULL.
}
\usage{getPrevAlign(object)}
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
loadIndex(reader,idx_file)
index.initialized(reader)
coords<-as.integer(c(0,0,1874))
range<-bamRange(reader,coords)
align<-getPrevAlign(range)
bamClose(reader)
}
\keyword{bamAlign}
\keyword{bamRange}
