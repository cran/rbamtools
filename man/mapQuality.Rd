\name{mapQuality}
\alias{mapQuality}
\title{mapQuality: Reads the mapQuality value from bamAlign.}
\description{The mapQuality is given by the alignment algorithm which produces the alignments.}
\usage{mapQuality(object)}
\arguments{
  \item{object}{An instance of class bamAlign.}
}
\details{mapQuality gives measures for the probability that the map-Position is wrong. A value 255 indicates that the mapping quality is not available.}
\value{Integer}
\references{
\verb{The SAM Format Specification. The SAM Format Specification Working Group. September 7, 2011.}
\url{http://www.samtools.sourceforge.net/SAM1.pdf}
}
\author{Wolfgang Kaisers}
\examples{
bam_file <- system.file("extdata", "accepted_hits.bam", package="rbamtools")
reader<-bamReader(bam_file)
isOpen(reader)
align<-getNextAlign(reader)
mapQuality(align)
bamClose(reader)
}
\keyword{bamAlign}
\keyword{mapQuality}
