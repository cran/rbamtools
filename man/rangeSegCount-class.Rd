\name{rangeSegCount-class}
\Rdversion{1.1}
\docType{class}
\alias{rangeSegCount-class}
%
\alias{rangeSegCount}
%
\alias{as.data.frame.rangeSegCount}
\alias{coerce, rangeSegCount, data.frame-method}
\alias{initialize, rangeSegCount-method}
%
\alias{meltDownSegments}
\alias{meltDownSegments-methods}
\alias{meltDownSegments,rangeSegCount-method}
%
\alias{size,rangeSegCount-method}
%
\title{Class \code{"rangeSegCount"}: Represents alignment counts in 
segments of genetic regions.}
%
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
% Description
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
\description{rangeSegCount.}
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
% Creation of objects
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
\section{Objects from the Class}{
    Objects can be created by calls of the form 
    \code{range<-rangeSegCount(object, coords, segments, complex)}.
}
%
\section{Slots}{
    \describe{
        \item{\code{position}:}{Integer. Vector of chromosomal positions
            which define the segments in which alignments are counted.}
        %
        \item{\code{count}:}{Integer. Number of alignments which are located
            in the position segments}
        %
        \item{\code{refname}:}{Character. Name of the reference sequence
            as given in the header section of the BAM file (e.g. 'chr1').}
        %
        \item{\code{LN}:}{Integer. Total size of reference sequence
            (chromosome).}
        %
        \item{\code{coords}:}{Numeric. A numeric vector giving the
            genomic coordinates from which alignments were counted. }
        %
        \item{\code{complex}:}{Logical.}
    }
}
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
% Details
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
\details{Two kinds of coordinates are used: A) coords (numeric) containing
seqid, start and end position (see also bamRange) which are used for reading
alignments from a BAM file, and B) position (numeric) which contains an 
arbitrary sized vector. The position values define the borders of genomic
segments in which alingments are counted. The segments usually would
cover a whole chromosome and separates genetic and intergenic regions.
The result is then the number of alignments in genetic and intergenetic
regions which then can be used as gene expression estimates.
The segments are defined as right open intervals.
When the position is given by c(1,10,20) for example, the first segment
is defined be nucleotide positions 1 to 9 and the second segment is 10 to 19.
The counting mechanism only takes the alignment position into account.
A possible overhang of the alingment over a segment boundary is not
represented here (because this would complicate the model of 'counting'
whole alignments.)
}
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
% Methods
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
\section{Methods}{
    \describe{
        \item{as.data.frame}{\code{signature(x="rangeSegCount")}:
            Returns data.frame representation of aligns.}
        %
        \item{coerce}{\code{signature(from="rangeSegCount", to="data.frame")}: 
            Coercion of rangeSegCount to data.frame.}
        %
        \item{meltDownSegments}{\code{signature(object="rangeSegCount",
            factor="numeric")}: Defines the factor by which the
            number of counting segments is shrinked. factor=2 means
            that counts for two adjacent segments are accumulated.}
    }
}
\author{Wolfgang Kaisers}
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
% Examples
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
\examples{
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# A) Open reader
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
bam <- system.file("extdata", "accepted_hits.bam", package="rbamtools")
idx<-paste(bam,"bai",sep=".")
# Open BAM file
reader<-bamReader(bam, idx=TRUE)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# B) Count range segment
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
coords <- c(0, 0, 2e4)
segments <- seq(14000, 20000, 20)
segcount<-rangeSegCount(reader, coords, segments)
segcount
dfr<-as.data.frame(segcount)
\donttest{
plot(count~position, dfr, type="l", las=1, bty="n", lwd=1.5)
}
}
\keyword{classes}
\keyword{bamReader}
