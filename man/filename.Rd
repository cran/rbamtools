\name{filename}
\alias{filename}
\title{
filename(bamReader)
}
\description{Returns filename of opened BAM-file from bamReader object.}
\usage{filename(object)}
\arguments{
  \item{object}{An object of class bamReader.}
}
\value{Character string representing filename of BAM-file.}
\author{Wolfgang Kaisers}
\examples{
library(rbamtools)
bam_file <- system.file("extdata", "accepted_hits.bam", package="rbamtools")
reader<-bamReader(bam_file)
filename(reader)
bamClose(reader)
}
\keyword{bamReader}
\keyword{filename}
