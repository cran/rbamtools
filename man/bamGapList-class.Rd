\name{bamGapList-class}
\Rdversion{1.1}
\docType{class}
\alias{bamGapList-class}
\alias{as.data.frame.bamGapList}
\alias{summary.bamGapList}
\alias{merge.bamGapList}
\alias{coerce,bamGapList,data.frame-method}
\alias{bamGapList}
\alias{bamGapList-methods}
\alias{bamGapList,bamReader-method}
\alias{initialize,bamGapList-method}
\alias{nAligns,bamGapList-method}
\alias{nAlignGaps,bamGapList-method}
\alias{show,bamGapList-method}
\alias{size,bamGapList-method}
\title{Class \code{"bamGapList"}}
\description{
The bamGapList class represents a list of Alignment gap (N-items in Cigar-data) sites . For each gap-site, left and right start and end positions as well as the gap-length are reported. Numbers of aligns supporting this site, number of left-sided start positions (<= 8) and the sum of overlapping nucleotides on the left side are given.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{siteList(reader,coords))}.
}
\section{Slots}{
  \describe{
    \item{\code{list}:}{\code{"externalptr"}. Point to double linked list struct.}
    \item{\code{refdata}:}{\code{"data.frame"}. Contains bamHeader like data for stored aligns.}
  }
}
\section{Methods}{
  \describe{
    \item{size}{\code{signature(x = "bamGapList")}: Returns number of site-items in list. }
    \item{coerce}{\code{signature(from = "bamGapList", to = "data.frame")}: Coercion of bamGapList to data.frame.}
    \item{coerce}{\code{signature("bamGapList","data.frame")}: Coercion of bamGapList to data.frame.}
    \item{show}{\code{signature(object = "bamGapList")}: Prints a short message with some summarizing data. }
    \item{nAligns}{\code{signature(object = "bamGapList")}: Returns number of aligns in specified Range. }
    \item{nAlignGaps}{\code{signature(object = "bamGapList")}: Returns number of align gaps in specified Range. }
	 }
}
\author{Wolfgang Kaisers}
\examples{
bam<-system.file("extdata", "accepted_hits.bam", package="rbamtools")
reader<-bamReader(bam,idx=TRUE)
bsl<-bamGapList(reader)
bsl
size(bsl)
nAligns(bsl)
nAlignGaps(bsl)
summary(bsl)
dfr<-as.data.frame(bsl)
head(dfr)
bamClose(reader)
}
\keyword{bamGapList}
