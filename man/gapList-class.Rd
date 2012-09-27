\name{gapList-class}
\Rdversion{1.1}
\docType{class}
\alias{as.data.frame,gapList-method}
\alias{gapList-class}
\alias{gapList,bamReader-method}
\alias{gapList}
\alias{initialize,gapList-method}
\alias{show,gapList-method}
\alias{size,gapList-method}
\title{Class \code{"gapList"}}
\description{
The gapList class represents a list of Alignment gaps (i.e. N-items in Cigar-data). For each gap, the type of left and right adjacent Cigar items is reported plus size of each. The list can be converted into a data.frame which then contains the columns: refid, position (which identify the align), left_cigar_len,left_cigar_type, left_stop (characterize the left boundary) and right_start, right_cigar_len, right_cigar_type (characterize the right boundary). The adjacent cigar-types should be 0 (i.e. M=match).
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("bamAlign",reader,coords))}.
}
\section{Slots}{
  \describe{
    \item{\code{list}:}{\code{"externalptr"}. Point to double linked list struct.}
  }
}
\section{Methods}{
  \describe{
    \item{size}{\code{signature(object = "gapList")}: Returns number of gapped-align items in list. }
    \item{as.data.frame}{\code{signature(object = "gapList")}: Returns data contained in this list as data.frame. }
    \item{show}{\code{signature(object = "gapList")}: Returns head of data.frame representation. }
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
loadIndex(reader,idx)
coords<-as.integer(c(0,0,249250621))
gl<-gapList(reader,coords)
df<-as.data.frame(gl)
head(df)
bamClose(reader)
}
\keyword{bamAlign}
