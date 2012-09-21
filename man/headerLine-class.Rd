\name{headerLine-class}
\Rdversion{1.1}
\docType{class}
\alias{headerLine-class}
\alias{initialize,headerLine-method}
\alias{getHeaderText,headerLine-method}
\alias{getVal,headerLine-method}
\alias{setVal,headerLine-method}
\alias{as.list,headerLine-method}
\title{Class \code{"headerLine"}}
\description{headerLine represents the header Line segment of header section.}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{hl<-new("headerLine")}.}
\section{Slots}{
  \describe{
    \item{\code{SO}:}{character}
    \item{\code{VN}:}{character}
  }
}
\section{Methods}{
  \describe{
    \item{getHeaderText}{\code{signature(x = "bamRange")}: ... }
    \item{getVal}{\code{signature(object = "bamRange")}: ... }
    \item{setVal}{\code{signature(object = "bamRange")}: ... }
    \item{as.list}{\code{signature(.Object = "bamRange")}: ... }
	 }
}
\author{Wolfgang Kaisers}
\examples{
bam<-system.file("extdata","accepted_hits.bam",package="rbamtools")
reader<-bamReader(bam)
header<-getHeader(reader)
htxt<-getHeaderText(header)
headLine<-headerLine(htxt)
getVal(headLine,"VN")
getVal(headLine,"SO")
setVal(headLine,"SO","unsorted")
l<-as.list(headLine)
txt<-getHeaderText(headLine)
bamClose(reader)
}
\keyword{classes}
\keyword{bamReader}
\keyword{bamHeader}
\keyword{headerLine}