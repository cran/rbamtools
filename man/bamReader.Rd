\name{bamReader}
\alias{bamReader}
\title{bamReader}
\description{The bamReader function returns an object of class bamReader. The bamReader object represents a connection to an opened BAM-file.}
\usage{bamReader(filename)}
\arguments{
  \item{filename}{Filename of BAM-file to be opened for reading.}
}
\author{Wolfgang Kaisers}
\examples{
bam<-system.file("extdata","accepted_hits.bam",package="rbamtools")
reader<-bamReader(bam)
align<-getNextAlign(reader)
name(align)
position(align)
bamClose(reader)
}
\keyword{bamReader}
