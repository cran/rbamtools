\name{matePosition}
\alias{matePosition}
\alias{mateRefID}
\title{Read position and reference-ID of mate alignment in paired end reads from a bamAlign.}
\description{Reads the position and reference-ID of mate alignment in paired end reads from a bamAlign.}
\usage{
matePosition(object)
mateRefID(object)
}
\arguments{
  \item{object}{An instance of class bamAlign.}
}
\value{
matePosition, mateRefID: Integer.
mateReverseStand, mateUnmapped: Logical.
}
\references{
\verb{The SAM Format Specification (v1.4-r985). The SAM Format Specification Working Group. September 7, 2011.}
\url{http://www.samtools.sourceforge.net/SAM1.pdf}
}
\author{Wolfgang Kaisers}
\examples{
bam_file <- system.file("extdata", "accepted_hits.bam", package="rbamtools")
reader<-bamReader(bam_file)
isOpen(reader)
align<-getNextAlign(reader)
matePosition(align)
mateRefID(align)
bamClose(reader)
}
\keyword{matePosition}
\keyword{mateRefID}
\keyword{bamAlign}
