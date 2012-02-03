\name{getRefData}
\alias{getRefData}
\title{getRefData: Returns an unparsed character string containing the header information of a BAM-file.}
\description{The header section contains various data fields which provides metadata about the stored alignment results.}
\usage{getRefData(object)}
\arguments{
  \item{object}{An object of class bamReader.}
}
\value{
data.frame
The returned data.frame contains three columns: For each reference sequence, the corresponding row contains the Reference-ID (1st column, refID), the Reference name (2nd column, refName) and the length of the Reference sequence (3rd column, refLength).
}
\references{
\verb{The SAM Format Specification (v1.4-r985). The SAM Format Specification Working Group. September 7, 2011.}
\url{http://www.samtools.sourceforge.net/SAM1.pdf}
}
\author{Wolfgang Kaisers}
\seealso{getHeaderText}
\examples{
bam_file <- system.file("extdata", "accepted_hits.bam", package="rbamtools")
reader<-bamReader(bam_file)
if(isOpen(reader))
{
ref<-getRefData(reader)
}
bamClose(reader)
}
\keyword{bamReader}
\keyword{getRefData}
