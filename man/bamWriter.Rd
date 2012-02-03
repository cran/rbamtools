\name{bamWriter}
\alias{bamWriter}
\title{bamWriter Object}
\description{bamWriter represents a writing connection to a BAM-file}
\usage{
bamWriter(reader,filename)
}
\arguments{
  \item{reader}{bamWriter imports header data (Reference table and metadata) from reader.}
  \item{filename}{Filename of BAM-file to be opened for writing.}
}
\author{Wolfgang Kaisers}

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

writer<-bamWriter(reader,"test.bam")
bamSave(writer,align)
bamSave(writer,range)

bamClose(reader)
bamClose(writer)
}
\keyword{bamWriter}
\keyword{bamAlign}
