\name{bamClose}
\alias{bamClose}
\title{bamClose(bamReader), bamClose(bamWriter)}
\description{Closes open connection from bamReader (or bamWriter) to .bam file.}
\usage{bamClose(object)}
\arguments{
  \item{object}{An instance of class bamReader or bamWriter.}
}
\details{
Closing a connection also deletes the stored index.
A closed connection can only be reopened by creating a new bamReader (or bamWriter) object.
}
\author{Wolfgang Kaisers}
\examples{
bam_file <- system.file("extdata", "accepted_hits.bam", package="rbamtools")
reader<-bamReader(bam_file)
align<-getNextAlign(reader)
bamClose(reader)
}
\keyword{bamClose}
\keyword{bamOpen}
