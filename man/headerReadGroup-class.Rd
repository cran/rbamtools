\name{headerReadGroup-class}
\Rdversion{1.1}
\docType{class}
%
\alias{addReadGroup}
\alias{addReadGroup-methods}
\alias{addReadGroup,headerReadGroup-method}
%
\alias{headerReadGroup-class}
%
\alias{initialize,headerReadGroup-method}
%
\alias{as.list,headerReadGroup-method}
%
\alias{getHeaderText,headerReadGroup-method}
%
\alias{getVal,headerReadGroup-method}
%
\alias{setVal,headerReadGroup-method}
%
\title{Class \code{"headerReadGroup"}: Representation of read - group data 
    in header section of BAM files}
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
% Description
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
\description{headerReadGroup represents data contained in the header
    section of BAM-files.}
\section{Objects from the Class}{
    Objects can be created by calls of the form
        \code{hl<-new("headerReadGroup")}.}
\section{Slots}{
    \describe{
        \item{\code{nrg}:}{\code{integer}. Number of read groups.}
        %
        \item{\code{ntl}:}{\code{integer}. 
            Number of different read group tags (currently = 12). }
        %
        \item{\code{ID}:}{\code{character}. Read group identifier.
            Each read group must have an unique ID.}
        %
        \item{\code{CN}:}{\code{character}. Name of sequencing center.}
        %
        \item{\code{DS}:}{\code{character}. Description.}
        %
        \item{\code{DT}:}{\code{character}. Date of run.}
        %
        \item{\code{FO}:}{\code{character}. Flow order.}
        %
        \item{\code{KS}:}{\code{character}. Array of nucleotide bases.}
        %
        \item{\code{LB}:}{\code{character}. Library.}
        %
        \item{\code{PG}:}{\code{character}. Programs used for processing.}
        %
        \item{\code{PI}:}{\code{character}. Predicted median insert size.}
        %
        \item{\code{PL}:}{\code{character}. 
            Platform/technology for production of reads.}
        %
        \item{\code{PU}:}{\code{character}. Unique platform identifier.}
        %
        \item{\code{SM}:}{\code{character}. Sample name.}
    }
}
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
% Methods
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
\section{Methods}{
    \describe{
        \item{addReadGroup}{\code{signature(object="headerReadGroup",
            l="list")}: Adds new read group to object.}
        %
        \item{as.list}{\code{signature(x="headerReadGroup")}:
            Returns data stored in this object as list. }
        %
        \item{getVal}{\code{signature(object="headerReadGroup")}: 
            Returns value of given segment.}
        %
        \item{setVal}{\code{signature(object="headerReadGroup")}:
            Sets value of given segment.}
        %
        \item{getHeaderText}{\code{signature(.Object="headerReadGroup")}:
            Returns textual representation of data as specified in SAM Format.}
}}
\references{
\verb{The SAM Format Specification (v1.4-r985). The SAM Format Specification Working Group. September 7, 2011.}
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
readGroup<-headerProgram(htxt)
getVal(readGroup,"ID")
setVal(readGroup,"ID","newID")
l<-as.list(readGroup)
txt<-getHeaderText(readGroup)
bamClose(reader)
}
\keyword{classes}
\keyword{bamReader}
\keyword{headerReadGroup}
