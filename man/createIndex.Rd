\name{createIndex}
\alias{createIndex}
\title{ createIndex(bamReader,filename) }
\description{ Creates an BAM-Index-file from existing BAM-file.}
\usage{ createIndex(object,idx_filename) }
\arguments{
  \item{object}{ An Instance of class bamReader }
  \item{idx_filename}{ A filename for new BAM-Index-file (usually bam_filename.bam.idx)}
}
\details{The function extracts the filename of the BAM-file.}
\value{None }
\author{Wolfgang Kaisers}
\examples{
library(rbamtools)
bam_file <- system.file("extdata", "accepted_hits.bam", package="rbamtools")
reader<-bamReader(bam_file)
idx_file<-"index.bam.bai"
createIndex(reader,idx_file)
loadIndex(reader,idx_file)
index.initialized(reader)
bamClose(reader)
}
\keyword{bamReader}
\keyword{Index}
