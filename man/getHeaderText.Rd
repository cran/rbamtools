\name{getHeaderText}
\alias{getHeaderText}
\title{Reads Header section Text of BAM-file.}
\description{
The header section contains various data fields which provide metadata about the stored alignment results. Returns an unparsed character string containing the header section of a BAM-file. Part of the Header-Text is the RefData information.
}
\usage{getHeaderText(object, delim = "\t")}
\arguments{
  \item{object}{An object of class bamReader}
  \item{delim}{Character used as tag delimiter}
}
\details{
The Header section of a BAM-file consists of the header line and the Reference sequence dictionary.
The header line contains information about the format version and a specification of the sorting order of alignments.
The reference dictionary contains information about the name and the length of each reference sequence.
}
\value{Character string.}
\author{Wolfgang Kaisers}
\references{
\verb{The SAM Format Specification (v1.4-r985). The SAM Format Specification Working Group. September 7, 2011.}
\url{http://www.samtools.sourceforge.net/SAM1.pdf}
}
\seealso{getRefData}
\examples{
library(rbamtools)
bam_file <- system.file("extdata", "accepted_hits.bam", package="rbamtools")
idx_file <- system.file("extdata", "accepted_hits.bam.bai", package="rbamtools")
reader<-bamReader(bam_file)
if(isOpen(reader))
	header<-getHeaderText(reader)
# ...
bamClose(reader)
}
\keyword{bamReader}
\keyword{bamHeader}
