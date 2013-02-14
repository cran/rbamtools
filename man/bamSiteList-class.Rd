\name{bamSiteList-class}
\Rdversion{1.1}
\docType{class}
\alias{bamSiteList-class}
\alias{as.data.frame.bamSiteList}
\alias{summary.bamSiteList}
\alias{merge.bamSiteList}
\alias{coerce,bamSiteList,data.frame-method}
\alias{bamSiteList}
\alias{bamSiteList-methods}
\alias{bamSiteList,bamReader-method}
\alias{initialize,bamSiteList-method}
\alias{nAligns,bamSiteList-method}
\alias{nGapAligns,bamSiteList-method}
\alias{show,bamSiteList-method}
\alias{size,bamSiteList-method}
\title{Class \code{"bamSiteList"}}
\description{
The bamSiteList class represents a list of Alignment gap (N-items in Cigar-data) sites . For each gap-site, left and right start and end positions as well as the gap-length are reported. Numbers of aligns supporting this site, number of left-sided start positions (<= 8) and the sum of overlapping nucleotides on the left side are given.
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
    \item{size}{\code{signature(x = "bamSiteList")}: Returns number of site-items in list. }
    \item{coerce}{\code{signature(from = "bamSiteList", to = "data.frame")}: Coercion of bamSiteList to data.frame.}
    \item{coerce}{\code{signature("bamSiteList","data.frame")}: Coercion of bamSiteList to data.frame.}
    \item{show}{\code{signature(object = "bamSiteList")}: Prints a short message with some summarizing data. }
    \item{nAlign}{\code{signature(object = "bamSiteList")}: Returns number of aligns in specified Range. }
    \item{nGapAlign}{\code{signature(object = "bamSiteList")}: Returns number of gapped aligns in specified Range. }
	 }
}
\author{Wolfgang Kaisers}
\examples{
library(rbamtools)
bam<-system.file("extdata", "accepted_hits.bam", package="rbamtools")
reader<-bamReader(bam,idx=TRUE)
bsl<-bamSiteList(reader)
bsl
size(bsl)
nAligns(bsl)
nGapAligns(bsl)
summary(bsl)
dfr<-as.data.frame(bsl)
head(dfr)
bamClose(reader)
}
\keyword{bamSiteList}
