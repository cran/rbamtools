\name{headerProgram-class}
\Rdversion{1.1}
\docType{class}
\alias{headerProgram-class}
\alias{initialize,headerProgram-method}
%
\alias{as.list,headerProgram-method}
%
\alias{getVal,headerProgram-method}
\alias{setVal,headerProgram-method}
%
\alias{getHeaderText,headerProgram-method}
%
\title{Class \code{"headerProgram"}: Represention of header - program section
 of BAM header.}
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
% Description
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
\description{headerProgram represents data contained in the header section 
    of BAM-files.}
\section{Objects from the Class}{
    Objects can be created by calls of the form \code{hl<-new("headerProgram")}.}
\section{Slots}{
    \describe{
        \item{\code{l}:}{list}
    }
}
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
% Methods
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
\section{Methods}{
    \describe{
        \item{as.list}{\code{signature(x = "headerProgram")}: 
            Converts data in object into list.}
        %
        \item{getVal}{\code{signature(object = "headerProgram")}:
            Returns value of given Segment.}
        %
        \item{setVal}{\code{signature(object = "headerProgram")}: 
            Sets value of given segment.}
        %
        \item{getHeaderText}{\code{signature(.Object = "headerProgram")}:
            Returns textual representation of data as specified in SAM
            File Format.}
}}
\references{
\verb{The SAM Format Specification (v1.4-r985). 
The SAM Format Specification Working Group. September 7, 2011.}
\url{http://www.samtools.sourceforge.net/SAM1.pdf}
}
\author{Wolfgang Kaisers}
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
% Examples
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
\examples{
bam <- system.file("extdata", "accepted_hits.bam", package="rbamtools")
reader<-bamReader(bam)
isOpen(reader)
header<-getHeader(reader)
htxt<-getHeaderText(header)
headProg<-headerProgram(htxt)
getVal(headProg,"PN")
setVal(headProg,"PN","bwa")
l<-as.list(headProg)
ht<-getHeaderText(headProg)
bamClose(reader)
}
\keyword{classes}
\keyword{bamReader}
\keyword{headerProgram}
