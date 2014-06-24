\name{bamRange-class}
\Rdversion{1.1}
\docType{class}
\alias{bamRange-class}
\alias{as.data.frame.bamRange}
\alias{getNextAlign,bamRange-method}
\alias{coerce,bamRange,data.frame-method}
\alias{initialize,bamRange-method}
\alias{rewind,bamRange-method}
%
\alias{[-methods}
\alias{[,bamRange-method}
%
\alias{getAlignRange}
\alias{getAlignRange-methods}
\alias{getAlignRange,bamRange-method}
%
\alias{getCoords}
\alias{getCoords-methods}
\alias{getCoords,bamRange-method}
%
\alias{getParams}
\alias{getParams-methods}
\alias{getParams,bamRange-method}
%
\alias{getRefName}
\alias{getRefName-methods}
\alias{getRefName,bamRange-method}
%
\alias{getSeqLen}
\alias{getSeqLen-methods}
\alias{getSeqLen,bamRange-method}
%
\alias{getPrevAlign}
\alias{getPrevAlign-methods}
\alias{getPrevAlign,bamRange-method}
%
\alias{insertPastCurrent}
\alias{insertPastCurrent-methods}
\alias{insertPastCurrent,bamRange-method}
%
\alias{insertPreCurrent}
\alias{insertPreCurrent-methods}
\alias{insertPreCurrent,bamRange-method}
%
\alias{moveCurrentAlign}
\alias{moveCurrentAlign-methods}
\alias{moveCurrentAlign,bamRange-method}
%
\alias{pop_back}
\alias{pop_back-methods}
\alias{pop_back,bamRange-method}
%
\alias{pop_front}
\alias{pop_front-methods}
\alias{pop_front,bamRange-method}
%
\alias{push_back}
\alias{push_back-methods}
\alias{push_back,bamRange-method}
%
\alias{push_front}
\alias{push_front-methods}
\alias{push_front,bamRange-method}
%
\alias{size}
\alias{size-methods}
\alias{size,bamRange-method}
%
\alias{stepNextAlign}
\alias{stepNextAlign-methods}
\alias{stepNextAlign,bamRange-method}
%
\alias{stepPrevAlign}
\alias{stepPrevAlign-methods}
\alias{stepPrevAlign,bamRange-method}
%
\alias{writeCurrentAlign}
\alias{writeCurrentAlign-methods}
\alias{writeCurrentAlign,bamRange-method}
%
\title{Class \code{"bamRange"}}
%
\description{bamRange represents a double linked list of bamAlign objects which overlap with a defined region in a BAM-file. The bamRange-function retrieves all alignments in the depicted Range from BAM-File into a bamRange object. A bamRange object maintains a double-linked list of aligns. The list keeps a pointer to a current align structure for iteration purposes. Addidionally there are some summarizing values stored (which are displayed by \code{show}) which describe the range inside the reference from which the \code{bamRange} object was read (seqid,qrBegin,qrEnd,complex) and some statistis (size,qSeqMinLen,qSeqMaxLen).
}
%
\section{Objects from the Class}{
Objects can be created by calls of the form \code{range<-bamRange(reader, coords)}.}
\section{Slots}{
  \describe{
    \item{\code{range}:}{External pointer. Points to double linked list of bamAligns.}
  }
}
%
\details{\code{bamRange} objects internally keep the following values:
\tabular{rll}{
  1 \tab seqid      \tab 0-based index of reference sequence\cr
  2 \tab qrBegin    \tab 0-based left boundary of query region (query range begin)\cr
  3 \tab qrEnd      \tab 0-based right boundary of query region (query range end)\cr
  4 \tab complex    \tab 0= all aligns included, 1= only aligns with n_cigar > 1 included \cr
  5 \tab rSeqLen    \tab Length of reference sequence\cr
  6 \tab qSeqMinLen \tab Minimum of query sequence length (= read length)\cr
  7 \tab qSeqMaxLen \tab Maximum of query sequence length (= read length)
}
For the \code{bamRange} class exists a rudimentary subsetting ('[') operator. '[' allows only for indexes > 0 and <= size(x). Index values are sorted in ascending order before values are extracted.
}
\section{Methods}{
  \describe{
    \item{as.data.frame}{\code{signature(x = "bamRange")}: Returns data.frame representation of aligns. }
    \item{coerce}{\code{signature(from = "bamRange", to = "data.frame")}: Coercion of bamRange to data.frame.}
    \item{bamSave}{\code{signature(object = "bamRange")}: Saves aligns stored in this list to BAM-file via a bamWriter object. }
    \item{getAlignRange}{\code{signature(object = "bamRange")}: Iterates through the list and returns 0-based position of the leftmost and rightmost matching nucleotide in range.}
    \item{getNextAlign}{\code{signature(object = "bamRange")}: Returns next align from current position and shifts current position to next one. }
    \item{getParams}{\code{signature(object = "bamRange")}: Returns named vector of stored parameters}
    \item{getPrevAlign}{\code{signature(object = "bamRange")}: Returns previous align from current position and shifts current position to previous one. }
    \item{getRefName}{\code{signature(.Object = "bamRange")}: Returns the reference sequence name from which the range was retrieved. }
    \item{getQualDf}{\code{signature(object = "bamRange", prob="logical")}: Returns position dependent counts of phred quality values. }
    \item{getQualQuantiles}{\code{signature(object = "bamRange", quantiles = "numeric")}: Returns position dependent quantile values for phred scores.}
    \item{plotQualQuant}{\code{signature(object = "bamRange")}: Plots phred quality quantiles for sequence positions. }
    \item{initialize}{\code{signature(.Object = "bamRange")}: Initializes bamRange object. }
    \item{insertPastCurrent}{\code{signature(object = "bamRange")}: Inserts align past current position into list. }
    \item{insertPreCurrent}{\code{signature(object = "bamRange")}: Insert align before current position into list. }
    \item{pop_back}{\code{signature(object = "bamRange")}: Removes last align from list. }
    \item{pop_front}{\code{signature(object = "bamRange")}: Removes first align from list. }
    \item{push_back}{\code{signature(object = "bamRange")}: Adds align at the end of the list. }
    \item{push_front}{\code{signature(object = "bamRange")}: Adds align at the front of the list. }
    \item{rewind}{\code{signature(object = "bamRange")}: Shifts current align to position before first align. }
    \item{size}{\code{signature(object = "bamRange")}: Returns number of aligns in list. }
    \item{writeCurrentAlign}{\code{signature(object = "bamRange")}: Overwrites current align with given align. }
	 }
}
\author{Wolfgang Kaisers}
\examples{
## + + + + + + + + + + + + + + + + + + + + + + + + + + + ##
## A) Open reader
## + + + + + + + + + + + + + + + + + + + + + + + + + + + ##
bam <- system.file("extdata", "accepted_hits.bam", package="rbamtools")
idx<-paste(bam,"bai",sep=".")
# Open BAM file
reader<-bamReader(bam)
# create.index(reader,idx)
# Load BAM index file
load.index(reader,idx)
index.initialized(reader)   # Should return 'TRUE'

## + + + + + + + + + + + + + + + + + + + + + + + + + + + ##
## B) Read range
## + + + + + + + + + + + + + + + + + + + + + + + + + + + ##
#  Find appropriate refid (=ID)
#  Returns a data.frame with three columns:
#  ID=refid, SN=Sequence Name, LN=Sequence length

rdf<-getRefData(reader)
head(rdf)

## + + + + + + + + + + + + + + + + + + + + + + + + + + + ##
# The sequence length (LN) also determines valid
# range for start and stop coordinates
# Invalid refid-, start- or stop-coordinates will
# release an error.
# coords: refid=0, start=0, stop=249250621
## + + + + + + + + + + + + + + + + + + + + + + + + + + + ##
coords<-as.integer(c(0,0,249250621))
range<-bamRange(reader,coords)
size(range)
range<-bamRange(reader,coords,complex=TRUE)
dfr<-as.data.frame(range)
size(range)
align<-getNextAlign(range)
cigarData(align)
rewind(range)

# while(!is.null(align))
# {
#   print(position(align))
#   align<-getNextAlign(range)
# }
bamClose(reader)

## + + + + + + + + + + + + + + + + + + + + + + + + + + + ##
## C) Get print message and some other values
## + + + + + + + + + + + + + + + + + + + + + + + + + + + ##
range
getCoords(range)
getSeqLen(range)
getParams(range)
getRefName(range)
getAlignRange(range)

## + + + + + + + + + + + + + + + + + + + + + + + + + + + ##
## D) Rudimentary subsetting
## + + + + + + + + + + + + + + + + + + + + + + + + + + + ##
range[1:5]
}
\keyword{classes}
\keyword{bamReader}
\keyword{bamAlign}
