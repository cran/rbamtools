\name{bamHeader-class}
\Rdversion{1.1}
\docType{class}
\alias{as.character,bamHeader-method}
\alias{as.character}
\alias{bamHeader-class}
\alias{bamHeader}
\alias{show,bamHeader-method}
\alias{bamWriter,bamHeader-method}
\alias{getHeaderText,bamHeader-method}
\alias{initialize,bamHeader-method}
\title{Class \code{"bamHeader"}}
\description{bamHeader represents Data contained in the header Section of BAM-files.}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{header<-new("bamHeader")}.}
\section{Slots}{
  \describe{
    \item{\code{header}:}{Object of class \code{"externalptr"}. Points to samtools bam_header_t struct.}
  }
}
\section{Methods}{
  \describe{
    \item{getHeaderText}{\code{signature(x = "bamHeader")}: Returns textual representation of data stored in this class as desrcribed in SAM Format Specification.}
	 }
}
\references{
\verb{The SAM Format Specification (v1.4-r985). The SAM Format Specification Working Group. September 7, 2011.}
\url{http://www.samtools.sourceforge.net/SAM1.pdf}
}
\author{Wolfgang Kaisers}
\examples{
bam<-system.file("extdata","accepted_hits.bam",package="rbamtools")
# Open bam file
reader<-bamReader(bam)
# Extract binary header structure
header<-getHeader(reader)
header
# Extract textual representation
htxt<-getHeaderText(header)
# Extract header line segment
hl<-headerLine(htxt)
# Get header program
hp<-headerProgram(htxt)
hp
# Re-create binary header structure
header2<-bamHeader(htxt)
# Use created structure for new BAM-file
writer<-bamWriter(header2,"test.bam")
bamClose(reader)
bamClose(writer)
}
\keyword{classes}
\keyword{bamReader}
\keyword{bamHeader}
