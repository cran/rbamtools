\name{headerReadGroup-class}
\Rdversion{1.1}
\docType{class}
\alias{headerReadGroup-class}
\alias{initialize,headerReadGroup-method}
\alias{as.list,headerReadGroup-method}
\alias{getHeaderText,headerReadGroup-method}
\alias{getVal,headerReadGroup-method}
\alias{setVal,headerReadGroup-method}
\title{Class \code{"headerReadGroup"}}
\description{headerReadGroup represents Data contained in the header Section of BAM-files.}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{hl<-new("headerReadGroup")}.}
\section{Slots}{
  \describe{
    \item{\code{l}:}{list}
  }
}
\section{Methods}{
  \describe{
    \item{as.list}{\code{signature(x = "headerReadGroup")}: Returns data stored in this object as list. }
    \item{getVal}{\code{signature(object = "headerReadGroup")}: Returns value of given segment. }
    \item{setVal}{\code{signature(object = "headerReadGroup")}: Sets value of given segment. }
    \item{getHeaderText}{\code{signature(.Object = "headerReadGroup")}: Returns textual representation of data as specified in SAM Format.}
}}
\references{
\verb{The SAM Format Specification (v1.4-r985). The SAM Format Specification Working Group. September 7, 2011.}
\url{http://www.samtools.sourceforge.net/SAM1.pdf}
}
\author{Wolfgang Kaisers}
\examples{
bam <- system.file("extdata", "accepted_hits.bam", package="rbamtools")
reader<-bamReader(bam)
isOpen(reader)
header<-getHeader(reader)
htxt<-getHeaderText(header)
readGroup<-headerProgram(htxt)
getVal(readGroup,"ID")
setVal(readGroup,"ID","newID")
l<-as.list(readGroup)
txt<-getHeaderText(readGroup)
bamClose(reader)
}
\keyword{classes}
\keyword{bamReader}
\keyword{headerReadGroup}
