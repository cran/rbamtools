\name{refSeqDict-class}
\Rdversion{1.1}
\docType{class}
\alias{refSeqDict-class}
\alias{initialize,refSeqDict-method}
\alias{addSeq,refSeqDict-method}
\alias{dim,refSeqDict-method}
\alias{[,refSeqDict-method}
\alias{[<-,refSeqDict-method}
\alias{head,refSeqDict-method}
\alias{tail,refSeqDict-method}
\alias{removeRows,refSeqDict-method}
\alias{getHeaderText,refSeqDict-method}
\alias{as.data.frame,refSeqDict-method}
\alias{addSeq}
\alias{head}
\alias{tail}
\alias{removeRows}
\title{Class \code{"refSeqDict"}}
\description{refSeqDict represents Data contained in the header Section of BAM-files.}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{hl<-new("refSeqDict")}.}
\section{Slots}{
  \describe{
    \item{\code{df}:}{data.frame}
  }
}
\section{Methods}{
  \describe{
    \item{dim}{\code{signature(x = "refSeqDict")}: ... }
    \item{as.data.frame}{\code{signature(object = "refSeqDict")}: ... }
    \item{removeRows}{\code{signature(object = "refSeqDict")}: ... }
    \item{addSeq}{\code{signature(.Object = "refSeqDict")}: ... }
    \item{head}{\code{signature(.Object = "refSeqDict")}: ... }
    \item{tail}{\code{signature(.Object = "refSeqDict")}: ... }
    \item{getHeaderText}{\code{signature(.Object = "refSeqDict")}: ... }
	 }
}
\author{Wolfgang Kaisers}
\examples{
bam_file <- system.file("extdata", "accepted_hits.bam", package="rbamtools")
reader<-bamReader(bam_file)
isOpen(reader)
header<-getHeader(reader)
refSeqs<-getRefSeqDict(header)
dim(refSeqs)
ht<-getHeaderText(refSeqs)
df<-as.data.frame(refSeqs)
refSeqs[1,1]
refSeqs[1,1]<-"newSeqName"
addSeq(refSeqs,SN="nextSeqName",LN=1000)
bamClose(reader)
}
\keyword{classes}
\keyword{bamReader}
\keyword{refSeqDict}
