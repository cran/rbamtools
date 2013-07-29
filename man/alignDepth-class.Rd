\name{alignDepth-class}
\Rdversion{1.1}
\docType{class}
\alias{alignDepth-class}
\alias{alignDepth}
\alias{alignDepth-methods}
\alias{alignDepth,bamRange-method}
\alias{getDepth}
\alias{getDepth-methods}
\alias{getDepth,alignDepth-method}
\alias{getParams,alignDepth-method}
\alias{getPos}
\alias{getPos-methods}
\alias{getPos,alignDepth-method}
\alias{plotAlignDepth}
\alias{plotAlignDepth-methods}
\alias{plotAlignDepth,alignDepth-method}

\alias{show,alignDepth-method}
\title{Class \code{"alignDepth"}}
\description{
The alignDepth class represents 
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{alignDepth(object,gap))} on bamRange objects. From the \code{bamRange} object, the range is extracted and for each nucleotide position whithin this range the numbers of align matches are calculated. When \code{alignDepth} is called wich \code{gap=TRUE}, the function counts aligns solely for gap-adjacent match regions (cigar-op's).
}
\section{Slots}{
  \describe{
    \item{\code{depth}:}{\code{"integer"}. Align depth values.}
    \item{\code{pos}:}{\code{"integer"}. Corresponding (1-based) positions}
    \item{\code{params}:}{\code{"numeric"}. }
    \item{\code{refname}:}{\code{"character"}. Name of reference sequence from which \code{bamRange} was extracted.}
  }
}
\details{
The \code{params} vector contains the following values:
\tabular{rll}{
  1 \tab seqid      \tab 0-based index of reference sequence\cr
  2 \tab qrBegin    \tab 0-based left boundary of query region (query range begin)\cr
  3 \tab qrEnd      \tab 0-based right boundary of query region (query range end)\cr
  4 \tab complex    \tab 0= all aligns included, 1= only aligns with n_cigar > 1 included \cr
  5 \tab rSeqLen    \tab Length of reference sequence\cr
  6 \tab qSeqMinLen \tab Minimum of query sequence length (= read length)\cr
  7 \tab qSeqMaxLen \tab Maximum of query sequence length (= read length)
}
}
\section{Methods}{
  \describe{
  \item{show}{\code{signature(object = "alignDepth")}: Prints a short message with some summarizing data. }
  \item{plotAlignDepth}{\code{signature(object = "alignDepth","...")}: Plots align depth in a line-plot.}
  \item{getDepth}{\code{signature(object = "alignDepth")}: Returns numeric align depth values. }
  \item{getPos}{\code{signature(object = "alignDepth")}: Returns numeric position values for align depth. }
  \item{getParams}{\code{signature(object = "alignDepth")}: Returns numeric parameter values. }
}}
\author{Wolfgang Kaisers}
\examples{
library(rbamtools)
# Open (indexed) BAM file
bam<-system.file("extdata", "accepted_hits.bam", package="rbamtools")
reader<-bamReader(bam,idx=TRUE)
# WASH7P coordinates
start<-16500
end<-16700
coords<-as.integer(c(0,start,end))
range<-bamRange(reader,coords)
bamClose(reader)
ad<-alignDepth(range)
ad
}
\keyword{bamAlign}
