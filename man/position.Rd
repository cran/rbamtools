\name{position}
\alias{position}
\alias{name}
\alias{refID}
\alias{readBases}
\alias{qualities}
\alias{cigarData}
\title{bamAlign accessor functions.}
\description{
position: Position of bamAlign.
name: Textual identifier of read sequence as written by the sequencing machine.
refID: ID of reference sequence for align.
readBases: Base sequence of aligned read in input sequence file (e.g. fasta or fastq).
qualities: Sequence of quality indicators in input sequence file (e.g. fasta or fastq).
cigarData: Details of alignment as reported from alignment algorithm.
}
\usage{
position(object)
name(object)
refID(object)
readBases(object)
qualities(object)
cigarData(object)
}
\arguments{
  \item{object}{
An instance of class bamAlign.
}
}
\details{
position (POS): 1-based left boundary position of align mapping.
refID: ID of reference sequence. To get the reference name, the ID can be translated with a table which contains the refID's and the corresponding reference names. This data.frame can be obtained with the getRefData.bamReader function.
Cigar data typically looks like 101M or 45M3501N56M. The numbers denote the cigar-length. The characters denote the cigar-type.
qualities: The SAM format specification says that qualities are reported as "ASCII of Phred-scaled base quality+33".
}
\value{
position, refid: Integer.
name,readBases,qualities: Character string.
cigarData: data.frame
\item{comp1 }{Numeric vector with cigar-length values}
\item{comp2 }{Character vector with cigar-type values}
}
\references{
\verb{The SAM Format Specification (v1.4-r985). The SAM Format Specification Working Group. September 7, 2011.}
\url{http://www.samtools.sourceforge.net/SAM1.pdf}
}
\author{Wolfgang Kaisers}
\examples{
library(rbamtools)
bam_file <- system.file("extdata", "accepted_hits.bam", package="rbamtools")
reader<-bamReader(bam_file)
isOpen(reader)
align<-getNextAlign(reader)
position(align)
name(align)
refID(align)
readBases(align)
qualities(align)
cigarData(align)
bamClose(reader)
}
\keyword{bamAlign}
\keyword{position}
\keyword{refID}
\keyword{cigarData}
\keyword{readBases}
\keyword{qualities}
