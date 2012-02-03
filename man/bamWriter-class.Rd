\name{bamWriter-class}
\Rdversion{1.1}
\docType{class}
\alias{bamWriter-class}
\alias{initialize,bamWriter-method}

\title{Class \code{"bamWriter"}}
\description{
Encapsulates functionality to write bamAligns to BAM-file.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{writer<-bamWriter(reader,filename)}.
}
\section{Slots}{
  \describe{
    \item{\code{filename}:}{Object of class \code{"character"} ~~ }
    \item{\code{writer}:}{Object of class \code{"bamWriter"} ~~ }
  }
}
\section{Methods}{
  \describe{
    \item{bamClose}{\code{signature(object = "bamWriter")}: ... }
    \item{bamSave}{\code{signature(object = "bamWriter")}: ... }
    \item{filename}{\code{signature(object = "bamWriter")}: ... }
    \item{initialize}{\code{signature(.Object = "bamWriter")}: ... }
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
writer<-bamWriter(reader,"test.bam")
bamSave(writer,align)
bamSave(writer,range)
bamClose(reader)
bamClose(writer)
}
\keyword{classes}
