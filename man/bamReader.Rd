\name{bamReader}
\alias{bamReader}
\alias{loadIndex}
\alias{index.initialized}
\title{bamReader}
\description{bamReader represents a connection to a BAM-file.}
\usage{bamReader(filename)}
\arguments{
  \item{filename}{Filename of BAM-file to be opened for reading.}
}
\author{Wolfgang Kaisers}
\examples{
bam_file<-system.file("extdata","accepted_hits.bam", package="rbamtools")
idx_file<-system.file("extdata","accepted_hits.bam.bai",package="rbamtools")
reader<-bamReader(bam_file)
align<-getNextAlign(reader)
name(align)
position(align)
loadIndex(reader,idx_file)
index.initialized(reader)
bamClose(reader)
}
\keyword{bamReader}
\keyword{bamAlign}
