\name{bamReader-class}
\Rdversion{1.1}
\docType{class}
\alias{bamReader-class}
\alias{createIndex,bamReader-method}
\alias{getHeader,bamReader-method}
\alias{getHeaderText,bamReader-method}
\alias{getNextAlign,bamReader-method}
\alias{getRefCount,bamReader-method}
\alias{getRefData,bamReader-method}
\alias{bamSort,bamReader-method}
\alias{loadIndex,bamReader-method}
\alias{initialize,bamReader-method}
\alias{index.initialized,bamReader-method}

\title{Class \code{"bamReader"}}
\description{bamReader encapulates functionality for reading of bamAligns from a BAM-file.}
\section{Objects from the Class}{Objects can be created by calls of the form \code{reader<-bamReader(filename)}.}
\section{Slots}{
  \describe{
    \item{\code{filename}:}{Character denoting name of BAM-file.}
    \item{\code{reader}:}{External pointer to opened BAM-file.}
    \item{\code{index}:}{External pointer to BAM index (used for retrieving bamRange objects from BAM-file.}
    \item{\code{indexName}:}{Name of opened index file.}
  }
}
\section{Methods}{
  \describe{
    \item{bamClose}{\code{signature(object = "bamReader")}: ... }
    \item{createIndex}{\code{signature(object = "bamReader")}: ... }
    \item{filename}{\code{signature(object = "bamReader")}: ... }
    \item{getHeaderText}{\code{signature(object = "bamReader")}: ... }
    \item{getNextAlign}{\code{signature(object = "bamReader")}: ... }
    \item{getRefCount}{\code{signature(object = "bamReader")}: ... }
    \item{getRefData}{\code{signature(object = "bamReader")}: ... }
    \item{bamSort}{\code{signature(object = "bamReader")}: ... }
    \item{index<-}{\code{signature(.Object = "bamReader")}: ... }
    \item{initialize}{\code{signature(.Object = "bamReader")}: ... }
	 }
}
\references{
\verb{The SAM Format Specification (v1.4-r985). The SAM Format Specification Working Group. September 7, 2011.}
\url{http://www.samtools.sourceforge.net/SAM1.pdf}
}
\author{Wolfgang Kaisers}
\examples{
library(rbamtools)
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
