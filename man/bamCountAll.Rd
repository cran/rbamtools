\name{bamCountAll}
\alias{bamCountAll}
\alias{bamCountAll-methods}
\alias{bamCountAll,bamReader-method}
\title{bamCountAll}
\description{The bamCountAll function takes a bamReader object and returns a data.frame. The data.frame contains the counts for all contained reference sequences.
}
\usage{bamCountAll(object,verbose=FALSE)}
\arguments{
  \item{object}{An instance of bamReader. Must be opened and contain initialized index.}
	\item{verbose}{Logical. Determines the amount of textual output during runtime.}
}
\details{The method returns a data.frame. Each row contains count data for one reference sequence. Each column contains the counts for one CIGAR-OP type ("M", "I","N","S","H","P","=","X"). Columns with the total number of aligns, the refid (ID) and the length of the reference sequence (LN), as retrieved by 'getRefData' are added.
}
\examples{
library(rbamtools)
bam <- system.file("extdata", "accepted_hits.bam", package="rbamtools")
reader<-bamReader(bam,idx=TRUE)
count<-bamCountAll(reader,verbose=TRUE)
count
bamClose(reader)
}
\keyword{bamCount}
\keyword{getRefData}
