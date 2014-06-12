\name{bamCount}
\alias{bamCount}
\alias{bamCount-methods}
\alias{bamCount,bamReader-method}
\title{bamCount}
\description{The bamCount function takes a bamReader object, a set of reference coordinates and the 'complex' argument and returns an integer vector.}
\usage{bamCount(object,coords)}
\arguments{
	\item{object}{An instance of bamReader. Must be opened and contain initialized index.}
	\item{coords}{Integer vector of length 3: coords=c(refid,start,stop)}
}
\details{The method returns integer vector of length 10. Entries 1 to 9 contain the number of CIGAR-OP items found in the given range. The 10th entry is the total number of aligns in the range. The returned vector is named. The names of the vector are "M", "I","N","S","H","P","=","X" and nAligns. The first 9 names are the abbreviations for CIGAR-OP items which are defined in the SAM file format reference.}
\examples{
bam <- system.file("extdata", "accepted_hits.bam", package="rbamtools")
reader<-bamReader(bam,idx=TRUE)
coords<-getRefCoords(reader,"chr1")
count<-bamCount(reader,coords)
count
bamClose(reader)
}
\keyword{bamCount}
\keyword{bamCountAll}
\keyword{index}
