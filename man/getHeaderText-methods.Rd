\name{getHeaderText-methods}
\docType{methods}
\alias{getHeaderText}
\alias{getHeaderText-methods}
\title{ \code{getHeaderText} Assembling member data into Header-Text}
\description{
The header section contains various data fields which provide metadata 
about the stored alignment results.
Returns an unparsed character string containing the header section
of a BAM-file.
Part of the Header-Text is the RefData information.
}
\usage{getHeaderText(object, delim = "\n")}
\arguments{
    \item{object}{An object of class bamHeader or bamHeaderText}
    \item{delim}{Character used as tag delimiter}
}
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
% Methods
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
\section{Methods}{
    \describe{
        \item{\code{signature(object = "bamHeader")}}{
            An object of class bamHeader
}
}}
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
% Details
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
\details{
The Header section of a BAM-file consists of the header line and
the Reference sequence dictionary.
The header line contains information about the format version and a
specification of the sorting order of alignments.
The reference dictionary contains information about the name and the length 
of each reference sequence.
The bamHeader object can be used to initialize an object of class bamWriter.
The bamHeaderText can be used to inspect the content of the header.
}
\value{An object of class bamHeaderText (from argument bamHeader) 
    and a string (from argument bamHeaderText)}
\author{Wolfgang Kaisers}
\references{
\verb{The SAM Format Specification (v1.4-r985). 
The SAM Format Specification Working Group. September 7, 2011.}
\url{http://www.samtools.sourceforge.net/SAM1.pdf}
}
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
% Examples
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
\examples{
bam<-system.file("extdata","accepted_hits.bam",package="rbamtools")
idx<-paste(bam,"bai",sep=".")
# Open BAM file
reader<-bamReader(bam)
header<-getHeader(reader)
htxt<-getHeaderText(header)
getHeaderText(htxt)
bamClose(reader)
}
\keyword{methods}
\keyword{getHeaderText}
