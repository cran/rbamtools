\name{bamHeader-class}
\Rdversion{1.1}
\docType{class}
\alias{bamHeader-class}
\alias{initialize,bamHeader-method}
\alias{getHeader,bamHeader-method}
\alias{getHeaderLine,bamHeader-method}
\alias{getHeaderProgram,bamHeader-method}
\alias{getHeaderReadGroup,bamHeader-method}
\alias{getRefSeqDict,bamHeader-method}
\alias{getHeaderText,bamHeader-method}
\alias{getHeader}
\alias{getHeaderLine}
\alias{getHeaderProgram}
\alias{getHeaderReadGroup}
\alias{getRefSeqDict}
\title{Class \code{"bamHeader"}}
\description{bamHeader represents Data contained in the header Section of BAM-files.}
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
    \item{getHeaderLine}{\code{signature(x = "bamHeader")}: ... }
    \item{getRefSeqDict}{\code{signature(object = "bamHeader")}: ... }
    \item{getHeaderReadGroup}{\code{signature(object = "bamHeader")}: ... }
    \item{getHeaderProgram}{\code{signature(.Object = "bamHeader")}: ... }
    \item{getHeaderText}{\code{signature(.Object = "bamHeader")}: ... }
	 }
}
\author{Wolfgang Kaisers}
\examples{
bam_file <- system.file("extdata", "accepted_hits.bam", package="rbamtools")
reader<-bamReader(bam_file)
isOpen(reader)
header<-getHeader(reader)
headLine<-getHeaderLine(header)
readGroup<-getHeaderReadGroup(header)
headProg<-getHeaderProgram(header)
headTxt<-getHeaderText(header)
bamClose(reader)
}
\keyword{classes}
\keyword{bamReader}
\keyword{bamHeader}
