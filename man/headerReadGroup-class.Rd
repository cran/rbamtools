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
    \item{as.list}{\code{signature(x = "headerReadGroup")}: ... }
    \item{getVal}{\code{signature(object = "headerReadGroup")}: ... }
    \item{setVal}{\code{signature(object = "headerReadGroup")}: ... }
    \item{getHeaderText}{\code{signature(.Object = "headerReadGroup")}: ... }
	 }
}
\author{Wolfgang Kaisers}
\examples{
bam_file <- system.file("extdata", "accepted_hits.bam", package="rbamtools")
reader<-bamReader(bam_file)
isOpen(reader)
header<-getHeader(reader)
readGroup<-getHeaderProgram(header)
getVal(readGroup,"ID")
setVal(readGroup,"ID","newID")
l<-as.list(readGroup)
txt<-getHeaderText(readGroup)
bamClose(reader)
}
\keyword{classes}
\keyword{bamReader}
\keyword{headerReadGroup}
