\name{flag}
\alias{flag}
\alias{flag<-}
\alias{paired}
\alias{paired<-}
\alias{properPair}
\alias{properPair<-}
\alias{unmapped}
\alias{unmapped<-}
\alias{mateUnmapped}
\alias{mateUnmapped<-}
\alias{reverseStrand}
\alias{reverseStrand<-}
\alias{mateReverseStrand}
\alias{mateReverseStrand<-}
\alias{firstInPair}
\alias{firstInPair<-}
\alias{secondInPair}
\alias{secondInPair<-}
\alias{secondaryAlign}
\alias{secondaryAlign<-}
\alias{failedQC}
\alias{failedQC<-}
\alias{pcrORopt_duplicate}
\alias{pcrORopt_duplicate<-}
\title{
Retrieving or setting flagged values of bamAlign objects.
}
\description{
In the BAM-align objects, two (of the four) bytes of the core-flag integer contain 11 binary managed flags. The other two bytes contain the number of cigar items reported in the current align.
The described (reading and writing) accessor functions address these flags individually while the flag function reads and writes all four bytes (=32 bits) at once. Therefore the flag function should only be used for copying all values from one existing align to another.

The binary maintained flags are (in ascending numerical value of the flag):
paired
properPair
unmapped
mateUnmapped
reverseStrand
mateReverseStrand
firstInPair
secondInPair
secondaryAlign
failedQC
pcrORopt_duplicate
}

\usage{
flag(object)
paired(object)
properPair(object)
unmapped(object)
mateUnmapped(object)
reverseStrand(object)
mateReverseStrand(object)
firstInPair(object)
secondInPair(object)
secondaryAlign(object)
failedQC(object)
pcrORopt_duplicate(object)

}
\arguments{
  \item{object}{ An instance of class bamReader }
}

\details{
bamAlign flag is a 32-Bit integer in which 16 bits represent binary flag values and 16 bits give the number of cigar items.
There are accessor functions for all individual flags. This function can be used to read or set all values simulaneously.
The strand information is returned as "Reverse" because that's the way it is stored in the BAM-file.
}
\value{
flag: Integer.
All other reading accessors: Logical.
All writing accessors: None.
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

i<-flag(align)
flag(align)<-i
paired(align)
paired(align)<-FALSE
properPair(align)
properPair(align)<-FALSE
unmapped(align)
unmapped(align)<-FALSE
mateUnmapped(align)
mateUnmapped(align)<-FALSE
reverseStrand(align)
reverseStrand(align)<-FALSE
mateReverseStrand(align)
mateReverseStrand(align)<-FALSE
firstInPair(align)
firstInPair(align)<-FALSE
secondInPair(align)
secondInPair(align)<-FALSE
secondaryAlign(align)
secondaryAlign(align)<-FALSE
failedQC(align)
failedQC(align)<-FALSE
pcrORopt_duplicate(align)
pcrORopt_duplicate(align)<-FALSE

bamClose(reader)
}
\keyword{bamReader}
\keyword{flag}
\keyword{paired}
\keyword{properPair}
\keyword{unmapped}
\keyword{mateUnmapped}
\keyword{reverseStrand}
\keyword{mateReverseStrand}
\keyword{firstInPair}
\keyword{secondInPair}
\keyword{secondaryAlign}
\keyword{failedQC}
\keyword{pcrORopt_duplicate}
