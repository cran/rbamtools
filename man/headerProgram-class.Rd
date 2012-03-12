\name{headerProgram-class}
\Rdversion{1.1}
\docType{class}
\alias{headerProgram-class}
\alias{initialize,headerProgram-method}
\alias{as.list,headerProgram-method}
\alias{getVal,headerProgram-method}
\alias{setVal,headerProgram-method}
\alias{getHeaderText,headerProgram-method}
\title{Class \code{"headerProgram"}}
\description{headerProgram represents Data contained in the header Section of BAM-files.}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{hl<-new("headerProgram")}.}
\section{Slots}{
  \describe{
    \item{\code{l}:}{list}
  }
}
\section{Methods}{
  \describe{
    \item{as.list}{\code{signature(x = "headerProgram")}: ... }
    \item{getVal}{\code{signature(object = "headerProgram")}: ... }
    \item{setVal}{\code{signature(object = "headerProgram")}: ... }
    \item{getHeaderText}{\code{signature(.Object = "headerProgram")}: ... }
	 }
}
\author{Wolfgang Kaisers}
\examples{
bam_file <- system.file("extdata", "accepted_hits.bam", package="rbamtools")
reader<-bamReader(bam_file)
isOpen(reader)
header<-getHeader(reader)
headProg<-getHeaderProgram(header)
getVal(headProg,"PN")
setVal(headProg,"PN","bwa")
l<-as.list(headProg)
ht<-getHeaderText(headProg)
bamClose(reader)
}
\keyword{classes}
\keyword{bamReader}
\keyword{headerProgram}
