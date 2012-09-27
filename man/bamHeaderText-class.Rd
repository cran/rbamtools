\name{bamHeaderText-class}
\Rdversion{1.1}
\docType{class}
\alias{bamHeader,bamHeaderText-method}
\alias{bamHeaderText-class}
\alias{bamHeaderText}
\alias{getHeader,bamHeaderText-method}
\alias{getHeader}
\alias{getHeaderText,bamHeaderText-method}
\alias{headerLine,bamHeaderText-method}
\alias{headerLine}
\alias{headerLine<-,bamHeaderText-method}
\alias{headerLine<-}
\alias{headerProgram,bamHeaderText-method}
\alias{headerProgram}
\alias{headerProgram<-,bamHeaderText-method}
\alias{headerProgram<-}
\alias{headerReadGroup,bamHeaderText-method}
\alias{headerReadGroup}
\alias{headerReadGroup<-,bamHeaderText-method}
\alias{headerReadGroup<-}
\alias{initialize,bamHeaderText-method}
\alias{refSeqDict,bamHeaderText-method}
\alias{refSeqDict}
\alias{refSeqDict<-,bamHeaderText-method}
\alias{refSeqDict<-}
\title{Class \code{"bamHeader"}}
\description{bamHeader manages textual representation of data contained in the header section of BAM-files.The binary representation (from which new BAM-files can be created) is encapsulated in class 'bamHeader'. Both types can be converted into each other. }
\section{Objects from the Class}{
Objects can be created by calls of the form \code{hl<-new("bamHeader")}.}
\section{Slots}{
  \describe{
    \item{\code{com}:}{character}
    \item{\code{head}:}{headerLine}
    \item{\code{dict}:}{refSeqDict}
    \item{\code{group}:}{headerReadGroup}
    \item{\code{prog}:}{headerProgram}
  }
}
\section{Methods}{
  \describe{
    \item{headerLine}{\code{signature(x = "bamHeader")}: Gets headerLine object. }
    \item{headerLine<-}{\code{signature(x = "bamHeader")}: Sets headerLine object. }
    \item{refSeqDict}{\code{signature(object = "bamHeader")}: Gets refSeqDict object. }
    \item{refSeqDict<-}{\code{signature(object = "bamHeader")}: Sets refSeqDict object. }
    \item{headerReadGroup}{\code{signature(object = "bamHeader")}: Gets headerReadGroup object.  }
    \item{headerReadGroup<-}{\code{signature(object = "bamHeader")}: Sets headerReadGroup object. }
    \item{headerProgram}{\code{signature(.Object = "bamHeader")}: Gets headerProgram object. }
    \item{headerProgram<-}{\code{signature(.Object = "bamHeader")}: Sets headerProgram object. }
    \item{getHeaderText}{\code{signature(.Object = "bamHeader")}: Returns whole information encoded in a character string as described in SAM Format Specification.  }
	 }
}
\author{Wolfgang Kaisers}
\examples{
bam<-system.file("extdata","accepted_hits.bam",package="rbamtools")
reader<-bamReader(bam)
header<-getHeader(reader)
bamClose(reader)
htxt<-getHeaderText(header)
headLine<-headerLine(htxt)
headerLine(htxt)<-headLine
readGroup<-headerReadGroup(htxt)
headerReadGroup(htxt)<-readGroup
headProg<-headerProgram(htxt)
headerProgram(htxt)<-headProg
headTxt<-getHeaderText(htxt)
}
\keyword{classes}
\keyword{bamReader}
\keyword{bamHeader}
