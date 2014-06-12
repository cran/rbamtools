\name{bamAlign-class}
\Rdversion{1.1}
\docType{class}
\alias{bamAlign-class}
%
\alias{alignQual}
\alias{alignQual-methods}
\alias{alignQual,bamAlign-method}
%alignQualVal
\alias{alignQualVal}
\alias{alignQualVal-methods}
\alias{alignQualVal,bamAlign-method}
%
\alias{alignSeq}
\alias{alignSeq-methods}
\alias{alignSeq,bamAlign-method}
%
\alias{cigarData}
\alias{cigarData-methods}
\alias{cigarData,bamAlign-method}
%
\alias{failedQC}
\alias{failedQC-methods}
\alias{failedQC,bamAlign-method}
%
\alias{failedQC<-}
\alias{failedQC<--methods}
\alias{failedQC<-,bamAlign-method}
%
\alias{firstInPair}
\alias{firstInPair-methods}
\alias{firstInPair,bamAlign-method}
%
\alias{firstInPair<-}
\alias{firstInPair<--methods}
\alias{firstInPair<-,bamAlign-method}
%
\alias{flag}
\alias{flag-methods}
\alias{flag,bamAlign-method}
%
\alias{flag<-}
\alias{flag<--methods}
\alias{flag<-,bamAlign-method}
%
\alias{initialize,bamAlign-method}
%
\alias{insertSize}
\alias{insertSize-methods}
\alias{insertSize,bamAlign-method}
%
\alias{mapQuality}
\alias{mapQuality-methods}
\alias{mapQuality,bamAlign-method}
%
\alias{matePosition}
\alias{matePosition-methods}
\alias{matePosition,bamAlign-method}
%
\alias{mateRefID}
\alias{mateRefID-methods}
\alias{mateRefID,bamAlign-method}
%
\alias{mateReverseStrand}
\alias{mateReverseStrand-methods}
\alias{mateReverseStrand,bamAlign-method}
%
\alias{mateReverseStrand<-}
\alias{mateReverseStrand<--methods}
\alias{mateReverseStrand<-,bamAlign-method}
%
\alias{mateUnmapped}
\alias{mateUnmapped-methods}
\alias{mateUnmapped,bamAlign-method}
%
\alias{mateUnmapped<-}
\alias{mateUnmapped<--methods}
\alias{mateUnmapped<-,bamAlign-method}
\alias{name}
\alias{name-methods}
\alias{name,bamAlign-method}
\alias{nCigar}
\alias{nCigar-methods}
\alias{nCigar,bamAlign-method}
\alias{paired}
\alias{paired-methods}
\alias{paired,bamAlign-method}
\alias{paired<-}
\alias{paired<--methods}
\alias{paired<-,bamAlign-method}
\alias{pcrORopt_duplicate}
\alias{pcrORopt_duplicate-methods}
\alias{pcrORopt_duplicate,bamAlign-method}
\alias{pcrORopt_duplicate<-}
\alias{pcrORopt_duplicate<--methods}
\alias{pcrORopt_duplicate<-,bamAlign-method}
\alias{position}
\alias{position-methods}
\alias{position,bamAlign-method}
\alias{properPair}
\alias{properPair-methods}
\alias{properPair,bamAlign-method}
\alias{properPair<-}
\alias{properPair<--methods}
\alias{properPair<-,bamAlign-method}
\alias{refID}
\alias{refID-methods}
\alias{refID,bamAlign-method}
\alias{reverseStrand}
\alias{reverseStrand-methods}
\alias{reverseStrand,bamAlign-method}
\alias{reverseStrand<-}
\alias{reverseStrand<--methods}
\alias{reverseStrand<-,bamAlign-method}
\alias{secondaryAlign}
\alias{secondaryAlign-methods}
\alias{secondaryAlign,bamAlign-method}
\alias{secondaryAlign<-}
\alias{secondaryAlign<--methods}
\alias{secondaryAlign<-,bamAlign-method}
\alias{secondInPair}
\alias{secondInPair-methods}
\alias{secondInPair,bamAlign-method}
\alias{secondInPair<--methods}
\alias{secondInPair<-}
\alias{secondInPair<-,bamAlign-method}
\alias{unmapped}
\alias{unmapped-methods}
\alias{unmapped,bamAlign-method}
\alias{unmapped<-}
\alias{unmapped<--methods}
\alias{unmapped<-,bamAlign-method}
\title{Class \code{"bamAlign"}}
\description{
The bamAlign class represents the content of a single alignment stored in a BAM-file. An instance of this class can be created by reading an Alignment from an object of class bamReader or bamRange. Accessor functions allow reading and writing of object contents. The alignment can be rewritten to a BAM-file via a bamWriter object or stored in a bamRange object.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("bamAlign", alignment)}.
}
\section{Slots}{
  \describe{
    \item{\code{align}:}{Object of class \code{"externalptr"}}
  }
}
\section{Methods}{
  \describe{
    \item{alignQual}{\code{signature(object = "bamAlign")}: Returns quality scores of read (one for each query-base). }
    \item{alignSeq}{\code{signature(object = "bamAlign")}: Returns base sequence (of read). }
    \item{cigarData}{\code{signature(object = "bamAlign")}: Returns data.frame with listed cigar items.}
    \item{failedQC}{\code{signature(object = "bamAlign")}: Gets flag item. }
    \item{failedQC<-}{\code{signature(object = "bamAlign")}: Sets flag item. }
    \item{firstInPair}{\code{signature(object = "bamAlign")}: Gets flag item. }
    \item{firstInPair<-}{\code{signature(object = "bamAlign")}: Sets flag item. }
    \item{flag}{\code{signature(object = "bamAlign")}: Retrieves 16-bit flag values which contains information of all flags in binary form.}
    \item{flag<-}{\code{signature(object = "bamAlign")}: Sets 16-bit flag value (all flags at once). }
    \item{initialize}{\code{signature(.Object = "bamAlign")}: Initializes bamAlign object }
    \item{insertSize}{\code{signature(object = "bamAlign")}: Returns insertSize value. }
    \item{mapQuality}{\code{signature(object = "bamAlign")}: Returns mapQuality value. }
    \item{matePosition}{\code{signature(object = "bamAlign")}: Returns mate Position value. }
    \item{mateRefID}{\code{signature(object = "bamAlign")}: Returns mateRefID value. }
    \item{mateReverseStrand}{\code{signature(object = "bamAlign")}: Gets flag item. }
    \item{mateReverseStrand<-}{\code{signature(object = "bamAlign")}: Sets flag item.}
    \item{mateUnmapped}{\code{signature(object = "bamAlign")}: Gets flag item. }
    \item{mateUnmapped<-}{\code{signature(object = "bamAlign")}: Sets flag item. }
    \item{name}{\code{signature(object = "bamAlign")}: Returns align name (read Identifier). }
    \item{nCigar}{\code{signature(object = "bamAlign")}: Returns number of cigar items in align. }
    \item{paired}{\code{signature(object = "bamAlign")}: Gets flag item. } 
    \item{paired<-}{\code{signature(object = "bamAlign")}: Sets flag item. }
    \item{pcrORopt_duplicate}{\code{signature(object = "bamAlign")}: Gets flag item. }
    \item{pcrORopt_duplicate<-}{\code{signature(object = "bamAlign")}: Sets flag item. }
    \item{position}{\code{signature(object = "bamAlign")}: Returns align position (on Ref-Sequence). }
    \item{properPair}{\code{signature(object = "bamAlign")}: Gets flag item. }
    \item{properPair<-}{\code{signature(object = "bamAlign")}: Sets flag item. }
    \item{refID}{\code{signature(object = "bamAlign")}: Returns (0-based) ID of Reference Sequence (As indicated by ID column returned by getRefData). }
    \item{reverseStrand}{\code{signature(object = "bamAlign")}: Gets flag item. }
    \item{reverseStrand<-}{\code{signature(object = "bamAlign")}: Sets flag item. }
    \item{secondaryAlign}{\code{signature(object = "bamAlign")}: Gets flag item. }
    \item{secondaryAlign<-}{\code{signature(object = "bamAlign")}: Sets flag item. }
    \item{secondInPair}{\code{signature(object = "bamAlign")}: Gets flag item. }
    \item{secondInPair<-}{\code{signature(object = "bamAlign")}: Sets flag item. }
    \item{unmapped}{\code{signature(object = "bamAlign")}: Gets flag item. }
    \item{unmapped<-}{\code{signature(object = "bamAlign")}: Sets flag item. }
	 }
}
\references{
\verb{The SAM Format Specification (v1.4-r985). The SAM Format Specification Working Group. September 7, 2011.}
\url{http://www.samtools.sourceforge.net/SAM1.pdf}
}
\author{Wolfgang Kaisers}
\examples{
# retrieve align from file
bam<-system.file("extdata", "accepted_hits.bam", package="rbamtools")
reader<-bamReader(bam)
align<-getNextAlign(reader)
bamClose(reader)
# use some accessors
name(align)
refID(align)
position(align)
nCigar(align)
cigarData(align)
alignSeq(align)
alignQual(align)
alignQualVal(align)
reverseStrand(align)
reverseStrand(align)<-TRUE
reverseStrand(align)
}
\keyword{bamAlign}
