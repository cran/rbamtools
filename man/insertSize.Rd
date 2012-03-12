\name{insertSize}
\alias{insertSize}
\title{insertSize: Reads insert size from bamAlign}
\description{
Reads the core.isize field from an alignment. Insert size should represent the distance between mate pair reads (but this is not explicitly declared in the format documentation).
}
\usage{insertSize(object)}
\arguments{
  \item{object}{An instance of class bamAlign.}
}
\value{Integer.}
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
insertSize(align)
bamClose(reader)
}
\keyword{insertSize}
\keyword{bamAlign}
