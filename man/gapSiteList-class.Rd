\name{gapSiteList-class}
\Rdversion{1.1}
\docType{class}
\alias{gapSiteList-class}
\alias{as.data.frame.gapSiteList}
\alias{merge.gapSiteList}
\alias{coerce,gapSiteList,data.frame-method}
\alias{siteList,bamReader-method}
\alias{siteList}
\alias{initialize,gapSiteList-method}
\alias{nAligns,gapSiteList-method}
\alias{nGapAligns,gapSiteList-method}
\alias{show,gapSiteList-method}
\alias{size,gapSiteList-method}
\title{Class \code{"gapSiteList"}}
\description{
The gapList class represents a list of Alignment gap (N-items in Cigar-data) sites . For each gap-site, left and right start and end positions as well as the gap-length are reported. Numbers of aligns supporting this site, number of left-sided start positions (<= 8) and the sum of overlapping nucleotides on the left side are given.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{siteList(reader,coords))}.
}
\section{Slots}{
  \describe{
    \item{\code{list}:}{\code{"externalptr"}. Point to double linked list struct.}
  }
}
\section{Methods}{
  \describe{
    \item{size}{\code{signature(x = "gapList")}: Returns number of site-items in list. }
    \item{coerce}{\code{signature(from = "gapList", to = "data.frame")}: Coercion of gapSiteList to data.frame.}
    \item{coerce}{\code{signature("gapList","data.frame")}: Coercion of gapSiteList to data.frame.}
    \item{show}{\code{signature(object = "gapList")}: Prints a short message with some summarizing data. }
    \item{nAlign}{\code{signature(object = "gapList")}: Returns number of aligns in specified Range. }
    \item{nGapAlign}{\code{signature(object = "gapList")}: Returns number of gapped aligns in specified Range. }
	 }
}
\author{Wolfgang Kaisers}
\examples{
library(rbamtools)
bam<-system.file("extdata", "accepted_hits.bam", package="rbamtools")
idx<-paste(bam,"bai",sep=".")
# Open BAM file
reader<-bamReader(bam)
# Load BAM index file
load.index(reader,idx)
coords<-as.integer(c(0,0,249250621))
sl<-siteList(reader,coords)
size(sl)
nAligns(sl)
nGapAligns(sl)
sl
df<-as.data.frame(sl)
head(df)
bamClose(reader)
}
\keyword{bamAlign}