\name{isOpen}
\alias{isOpen}
\title{isOpen(bamReader), isOpen(bamWriter)}
\description{Checks for open connection from bamReader (or bamWriter) to .bam file.}
\usage{isOpen(con,rw="")}
\arguments{
  \item{con}{An instance of class bamReader or bamWriter.}
  \item{rw}{This argument is included solely for compatibility with the function template for connections in package base. It's not evaluated internally.}
}
\details{The function acutally checks, if the externalptr (reader or writer) slot has nil-value.}
\author{Wolfgang Kaisers}
\examples{
library(rbamtools)
bam_file <- system.file("extdata", "accepted_hits.bam", package="rbamtools")
reader<-bamReader(bam_file)
isOpen(reader)
align<-getNextAlign(reader)
bamClose(reader)
isOpen(reader)
}
\keyword{bamClose}
\keyword{bamOpen}
