\name{bamRange-class}
\Rdversion{1.1}
\docType{class}
\alias{bamRange-class}
\alias{as.data.frame,bamRange-method}
\alias{getNextAlign,bamRange-method}
\alias{getPrevAlign,bamRange-method}
\alias{initialize,bamRange-method}
\alias{insertPastCurrent,bamRange-method}
\alias{insertPreCurrent,bamRange-method}
\alias{pop_back,bamRange-method}
\alias{pop_front,bamRange-method}
\alias{push_back,bamRange-method}
\alias{push_front,bamRange-method}
\alias{size,bamRange-method}
\alias{writeCurrentAlign,bamRange-method}
\alias{windBack,bamRange-method}

\title{Class \code{"bamRange"}}
\description{bamRange represents a list of bamAlign objects which overlap with a defined region in a BAM-file.}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{range<-bamRange(reader, coords)}.}
\section{Slots}{
  \describe{
    \item{\code{range}:}{External pointer to double linked list of bamAligns.}
  }
}
\section{Methods}{
  \describe{
    \item{as.data.frame}{\code{signature(x = "bamRange")}: ... }
    \item{bamSave}{\code{signature(object = "bamRange")}: ... }
    \item{getNextAlign}{\code{signature(object = "bamRange")}: ... }
    \item{initialize}{\code{signature(.Object = "bamRange")}: ... }
    \item{insertPastCurrent}{\code{signature(object = "bamRange")}: ... }
    \item{insertPreCurrent}{\code{signature(object = "bamRange")}: ... }
    \item{pop_back}{\code{signature(object = "bamRange")}: ... }
    \item{pop_front}{\code{signature(object = "bamRange")}: ... }
    \item{push_back}{\code{signature(object = "bamRange")}: ... }
    \item{push_front}{\code{signature(object = "bamRange")}: ... }
    \item{size}{\code{signature(object = "bamRange")}: ... }
    \item{writeCurrentAlign}{\code{signature(object = "bamRange")}: ... }
	 }
}
\author{Wolfgang Kaisers}
\examples{
bam_file <- system.file("extdata", "accepted_hits.bam", package="rbamtools")
idx_file <- system.file("extdata", "accepted_hits.bam.bai", package="rbamtools")
reader<-bamReader(bam_file)
isOpen(reader)
align<-getNextAlign(reader)
loadIndex(reader,idx_file)
index.initialized(reader)
coords<-as.integer(c(0,0,1874))
range<-bamRange(reader,coords)
align<-getNextAlign(range)
bamClose(reader)
}
\keyword{classes}
\keyword{bamReader}
\keyword{bamAlign}
