\name{refSeqDict-class}
\Rdversion{1.1}
\docType{class}
\alias{as.data.frame.refSeqDict}
\alias{coerce,refSeqDict,data.frame-method}
\alias{[,refSeqDict-method}
\alias{[<-,refSeqDict-method}
\alias{addSeq,refSeqDict-method}
\alias{addSeq}
\alias{as.data.frame,refSeqDict-method}
\alias{dim,refSeqDict-method}
\alias{getHeaderText,refSeqDict-method}
\alias{head,refSeqDict-method}
\alias{head}
\alias{initialize,refSeqDict-method}
\alias{refSeqDict-class}
\alias{removeSeqs,refSeqDict-method}
\alias{removeSeqs}
\alias{tail,refSeqDict-method}
\alias{tail}
\title{Class \code{"refSeqDict"}: Repesentation of data from reference 
    sequence dictionary in BAM file header.}
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
% Description
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
\description{refSeqDict represents Data contained in the header section
    of BAM-files.}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{hl<-new("refSeqDict")}.}
\section{Slots}{
    \describe{
        \item{\code{SN}:}{character}
        \item{\code{LN}:}{numeric}
        \item{\code{AS}:}{character}
        \item{\code{M5}:}{numeric}
        \item{\code{SP}:}{character}
        \item{\code{UR}:}{character}
    }
}
\section{Methods}{
    \describe{
        \item{dim}{\code{signature(x="refSeqDict")}:
            Returns dimension of data: number of rows and 6 columns.}
        %
        \item{as.data.frame}{\code{signature(object="refSeqDict")}: 
            Combines data of Slots to data.frame.}
        %
        \item{coerce}{\code{signature(from="bamRange", to="data.frame")}: 
            Coercion of refSeqDice to data.frame.}
        %
        \item{removeSeqs}{\code{signature(object="refSeqDict",rows="numeric")}:
            Removes reference sequence (i.e. one row from data.frame) entry. }
        %
        \item{addSeq}{\code{signature(.Object="refSeqDict",SN,LN,AS,M5,SP,UR)}:
            Adds reference sequence (i.e. one row in data.frame) entry.}
        %
        \item{head}{\code{signature(.Object="refSeqDict",n)}:
            Returns head of data.frame representation.}
        %
        \item{tail}{\code{signature(.Object="refSeqDict",n)}: 
            Returns tail of data.frame representation.}
        %
        \item{getHeaderText}{\code{signature(.Object="refSeqDict")}:
            Returns textual representation of data stored in this class as
            desrcribed in SAM Format Specification.}
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
refSeqs<-refSeqDict(htxt)
dim(refSeqs)
ht<-getHeaderText(refSeqs)
df<-as.data.frame(refSeqs)
addSeq(refSeqs,SN="nextSeqName",LN=1000)
removeSeqs(refSeqs,3)
bamClose(reader)
}
\keyword{classes}
\keyword{bamReader}
\keyword{refSeqDict}
