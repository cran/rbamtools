\name{gapSiteList-class}
\Rdversion{1.1}
\docType{class}
%
\alias{gapSiteList-class}
\alias{initialize,gapSiteList-method}
%
\alias{as.data.frame.gapSiteList}
\alias{merge.gapSiteList}
\alias{coerce,gapSiteList,data.frame-method}
%
\alias{siteList}
\alias{siteList,missing,missing-method}
\alias{siteList,numeric,missing-method}
\alias{siteList,bamReader,numeric-method}
\alias{siteList,bamRange,missing-method}
%
\alias{nAligns,gapSiteList-method}
\alias{nAlignGaps,gapSiteList-method}
%
\alias{refID,gapSiteList-method}
%
\alias{refID<-}
\alias{refID<--methods}
\alias{refID<-,gapSiteList,numeric-method}
%
\alias{show,gapSiteList-method}
\alias{size,gapSiteList-method}
\title{Class \code{"gapSiteList"}: Representation of genomic alignment gap
    positions.}
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
% Description
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
\description{
The gapSiteList class represents pooled alignment-gap sites.
Gap-sites are characterized by unique refid, left-end and right-start positions
in each list.
Numbers of aligns supporting this site, number of left-sided start positions 
(<= 8) and the sum of overlapping nucleotides on the left side are given.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{siteList(reader,coords))}.
}
\section{Slots}{
    \describe{
        \item{\code{list}:}{\code{"externalptr"}. Point to double linked list struct.}
    }
}
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
% Methods
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
\section{Methods}{
    \describe{
        \item{size}{\code{signature(x = "gapSiteList")}:
            Returns number of site-items in list.}
        %
        \item{coerce}{\code{signature(from = "gapSiteList", to = "data.frame")}:
            Coercion of gapSiteList to data.frame.}
        %
        \item{coerce}{\code{signature("gapSiteList","data.frame")}: 
            Coercion of gapSiteList to data.frame.}
        %
        \item{show}{\code{signature(object = "gapSiteList")}:
            Prints a short message with some summarizing data.}
        %
        \item{nAligns}{\code{signature(object = "gapSiteList")}:
            Returns number of aligns in specified Range.}
        %
        \item{nAlignGaps}{\code{signature(object = "gapSiteList")}:
            Returns number of align gaps in specified Range.}
        %
        \item{refID}{\code{signature(object="gapSiteList")}: 
            Returns refID from which \code{gapSiteList} has been retrieved.}
        %
        \item{refID<-}{\code{signature(object="gap",
            value="numeric")}: Sets flag item.}
    }
}
\author{Wolfgang Kaisers}
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
% Examples
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
\examples{
# Open (indexed) BAM file
bam<-system.file("extdata", "accepted_hits.bam", package="rbamtools")
reader<-bamReader(bam,idx=TRUE)
coords<-getRefCoords(reader,"chr1")
sl<-siteList(reader,coords)
size(sl)
nAligns(sl)
nAlignGaps(sl)
sl
refID(sl)
df<-as.data.frame(sl)
head(df)
# Create from bamRange:
br <- bamRange(reader, coords)
sl <- siteList(br)
#
bamClose(reader)
# Create empty list
sl <- siteList()
}
\keyword{bamAlign}
