\name{headerLine-class}
\Rdversion{1.1}
\docType{class}
%
\alias{headerLine-class}
%
\alias{initialize,headerLine-method}
%
\alias{getHeaderText,headerLine-method}
%
\alias{getVal,headerLine-method}
\alias{setVal,headerLine-method}
\alias{as.list,headerLine-method}
%
\title{Class \code{"headerLine"}: Representation of header line segment of 
 header section for BAM files.}
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
% Description
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
\description{headerLine represents the header Line segment of header section.}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{hl<-new("headerLine")}.}
\section{Slots}{
    \describe{
        \item{\code{SO}:}{character}
        %
        \item{\code{VN}:}{character}
    }
}
\section{Methods}{
    \describe{
        \item{getHeaderText}{\code{signature(x = "bamRange")}:
            Retrieve textual representation of header}
        %
        \item{getVal}{\code{signature(object = "bamRange")}: 
            Retrieving values for given item names.}
        %
        \item{setVal}{\code{signature(object = "bamRange")}: 
            Setting values for given item names.}
        %
        \item{as.list}{\code{signature(.Object = "bamRange")}: Coercing of
            data into a list}
}}
\author{Wolfgang Kaisers}
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
% Examples
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
\examples{
bam<-system.file("extdata","accepted_hits.bam",package="rbamtools")
reader<-bamReader(bam)
header<-getHeader(reader)
htxt<-getHeaderText(header)
headLine<-headerLine(htxt)
getVal(headLine,"VN")
getVal(headLine,"SO")
setVal(headLine,"SO","unsorted")
l<-as.list(headLine)
txt<-getHeaderText(headLine)
bamClose(reader)
}
\keyword{classes}
\keyword{bamReader}
\keyword{bamHeader}
\keyword{headerLine}
