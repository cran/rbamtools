\name{gapList-class}
\Rdversion{1.1}
\docType{class}
\alias{gapList-class}
\alias{as.data.frame.gapList}
\alias{coerce,gapList,data.frame-method}
\alias{gapList,bamReader-method}
\alias{gapList}
\alias{initialize,gapList-method}
\alias{nAligns}
\alias{nAligns-methods}
\alias{nAligns,gapList-method}
\alias{nAlignGaps}
\alias{nAlignGaps-methods}
\alias{nAlignGaps,gapList-method}
\alias{show,gapList-method}
\alias{size,gapList-method}
\title{Class \code{"gapList"}}
\description{
The gapList class represents a list of Alignment gaps (i.e. N-items in Cigar-data). For each gap, the type of left and right adjacent Cigar items is reported plus size of each. The list can be converted into a data.frame which then contains the columns: refid, position (which identify the align), left_cigar_len,left_cigar_type, left_stop (characterize the left boundary) and right_start, right_cigar_len, right_cigar_type (characterize the right boundary). The adjacent cigar-types should be 0 (i.e. M=match). left_stop is the 0-based position of last exon nucleotide, right_start is the 0-based position of the first exon nucleotide.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{gapList(reader,coords))}.
}
\section{Slots}{
  \describe{
    \item{\code{list}:}{\code{"externalptr"}. Point to double linked list struct.}
  }
}
\section{Methods}{
  \describe{
    \item{size}{\code{signature(x = "gapList")}: Returns number of gapped-align items in list. }
    \item{coerce}{\code{signature(from = "gapList", to = "data.frame")}: Coercion of gapList to data.frame.}
    \item{as.data.frame}{\code{signature(x="gapList",row.names=NULL,optional=FALSE)}: Coercion of gapList to data.frame.}
    \item{show}{\code{signature(object = "gapList")}: Prints a short message with some summarizing data. }
    \item{nAligns}{\code{signature(object = "gapList")}: Returns number of aligns in specified Range. }
    \item{nAlignGaps}{\code{signature(object = "gapList")}: Returns number of align gaps in specified Range. }
	 }
}
\author{Wolfgang Kaisers}
\examples{
library(rbamtools)
# Open (indexed) BAM file
bam<-system.file("extdata", "accepted_hits.bam", package="rbamtools")
reader<-bamReader(bam,idx=TRUE)
coords<-getRefCoords(reader,"chr1")
gl<-gapList(reader,coords)
size(gl)
nAligns(gl)
nAlignGaps(gl)
gl
dfr<-as.data.frame(gl)
# coerce
dfr<-as(gl,"data.frame")
head(df)
bamClose(reader)
}
\keyword{bamAlign}
