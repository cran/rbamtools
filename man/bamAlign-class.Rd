\name{bamAlign-class}
\Rdversion{1.1}
\docType{class}
\alias{bamAlign-class}
\alias{cigarData,bamAlign-method}
\alias{failedQC<-,bamAlign-method}
\alias{failedQC,bamAlign-method}
\alias{firstInPair<-,bamAlign-method}
\alias{firstInPair,bamAlign-method}
\alias{flag<-,bamAlign-method}
\alias{flag,bamAlign-method}
\alias{initialize,bamAlign-method}
\alias{insertSize,bamAlign-method}
\alias{mapQuality,bamAlign-method}
\alias{matePosition,bamAlign-method}
\alias{mateRefID,bamAlign-method}
\alias{mateReverseStrand<-,bamAlign-method}
\alias{mateReverseStrand,bamAlign-method}
\alias{mateUnmapped<-,bamAlign-method}
\alias{mateUnmapped,bamAlign-method}
\alias{name,bamAlign-method}
\alias{paired<-,bamAlign-method}
\alias{paired,bamAlign-method}
\alias{pcrORopt_duplicate<-,bamAlign-method}
\alias{pcrORopt_duplicate,bamAlign-method}
\alias{position,bamAlign-method}
\alias{properPair<-,bamAlign-method}
\alias{properPair,bamAlign-method}
\alias{qualities,bamAlign-method}
\alias{readBases,bamAlign-method}
\alias{refID,bamAlign-method}
\alias{reverseStrand<-,bamAlign-method}
\alias{reverseStrand,bamAlign-method}
\alias{secondaryAlign<-,bamAlign-method}
\alias{secondaryAlign,bamAlign-method}
\alias{secondInPair<-,bamAlign-method}
\alias{secondInPair,bamAlign-method}
\alias{unmapped<-,bamAlign-method}
\alias{unmapped,bamAlign-method}

\title{Class \code{"bamAlign"}}
\description{
The bamAlign class represents the content of a single alignment stored in a BAM-file. An instance of this class can be created by reading an Alignment from a BAM-file. Accessor functions allow reading and writing of object contents. The alignment can be rewritten to a BAM-file via a bamWriter object.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("bamAlign", alignment)}.
}
\section{Slots}{
  \describe{
    \item{\code{align}:}{Object of class \code{"externalptr"} ~~ }
  }
}
\section{Methods}{
  \describe{
    \item{cigarData}{\code{signature(object = "bamAlign")}: ... }
    \item{failedQC<-}{\code{signature(object = "bamAlign")}: ... }
    \item{failedQC}{\code{signature(object = "bamAlign")}: ... }
    \item{firstInPair<-}{\code{signature(object = "bamAlign")}: ... }
    \item{firstInPair}{\code{signature(object = "bamAlign")}: ... }
    \item{flag<-}{\code{signature(object = "bamAlign")}: ... }
    \item{flag}{\code{signature(object = "bamAlign")}: ... }
    \item{initialize}{\code{signature(.Object = "bamAlign")}: ... }
    \item{insertSize}{\code{signature(object = "bamAlign")}: ... }
    \item{mapQuality}{\code{signature(object = "bamAlign")}: ... }
    \item{matePosition}{\code{signature(object = "bamAlign")}: ... }
    \item{mateRefID}{\code{signature(object = "bamAlign")}: ... }
    \item{mateReverseStrand<-}{\code{signature(object = "bamAlign")}: ... }
    \item{mateReverseStrand}{\code{signature(object = "bamAlign")}: ... }
    \item{mateUnmapped<-}{\code{signature(object = "bamAlign")}: ... }
    \item{mateUnmapped}{\code{signature(object = "bamAlign")}: ... }
    \item{name}{\code{signature(object = "bamAlign")}: ... }
    \item{paired<-}{\code{signature(object = "bamAlign")}: ... }
    \item{paired}{\code{signature(object = "bamAlign")}: ... }
    \item{pcrORopt_duplicate<-}{\code{signature(object = "bamAlign")}: ... }
    \item{pcrORopt_duplicate}{\code{signature(object = "bamAlign")}: ... }
    \item{position}{\code{signature(object = "bamAlign")}: ... }
    \item{properPair<-}{\code{signature(object = "bamAlign")}: ... }
    \item{properPair}{\code{signature(object = "bamAlign")}: ... }
    \item{qualities}{\code{signature(object = "bamAlign")}: ... }
    \item{queryBases}{\code{signature(object = "bamAlign")}: ... }
    \item{refID}{\code{signature(object = "bamAlign")}: ... }
    \item{reverseStrand<-}{\code{signature(object = "bamAlign")}: ... }
    \item{reverseStrand}{\code{signature(object = "bamAlign")}: ... }
    \item{secondaryAlign<-}{\code{signature(object = "bamAlign")}: ... }
    \item{secondaryAlign}{\code{signature(object = "bamAlign")}: ... }
    \item{secondInPair<-}{\code{signature(object = "bamAlign")}: ... }
    \item{secondInPair}{\code{signature(object = "bamAlign")}: ... }
    \item{unmapped<-}{\code{signature(object = "bamAlign")}: ... }
    \item{unmapped}{\code{signature(object = "bamAlign")}: ... }
	 }
}
\author{Wolfgang Kaisers}
\examples{
library(rbamtools)
bam_file <- system.file("extdata", "accepted_hits.bam", package="rbamtools")
reader<-bamReader(bam_file)
align<-getNextAlign(reader)
name(align)
position(align)
bamClose(reader)
}
\keyword{bamAlign}
