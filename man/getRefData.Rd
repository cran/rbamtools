\name{getRefData}
%
\alias{getRefCoords}
\alias{getRefCoords-methods}
\alias{getRefCoords,bamReader,character-method}
\alias{getRefCoords,bamReader,factor-method}
%
\alias{getRefCount}
\alias{getRefCount-methods}
\alias{getRefCount,bamReader-method}
%
\alias{getRefData}
\alias{getRefData-methods}
\alias{getRefData,bamReader-method}
\alias{getRefData,bamWriter-method}
%
\alias{getRefId}
\alias{getRefId-methods}
\alias{getRefId,bamReader,character-method}
\alias{getRefId,bamWriter,character-method}
\alias{getRefId,bamReader,factor-method}
\alias{getRefId,bamWriter,factor-method}
%
\title{Retrieve reference sequence from a BAM file as data.frame}
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
% Description
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
\description{The four functions: \code{getRefCoords}, \code{getRefCount},
\code{getRefData} and \code{getRefId} provide reading access to data about
the present reference sequences in bamReader or bamWriter.}
\usage{getRefData(object)}
\arguments{
        \item{object}{\code{bamReader}.
            The reader must be opened (otherwise an error is thrown).}
}
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
% Details
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
\details{\code{getRefData} returns a \code{data.frame} with three columns
(ID, SN, LN).
ID is the (0-based index which must be given when a \code{bamRange} is
extracted.
SN is the name of the sequence (e.g. chr1 for UCSC). 
LN is the length of the reference sequence.
\code{getRefCount} returns the number of reference sequences. 
\code{getRefCoords} returns a vector of coordinates which can be used to extract all stored aligns for this sequence from the \code{bamReader} into a \code{bamRange} object.
}
\author{Wolfgang Kaisers}
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
% Examples
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
\examples{
bam<-system.file("extdata","accepted_hits.bam",package="rbamtools")
reader<-bamReader(bam,idx=TRUE)
getRefData(reader)
getRefCount(reader)
getRefId(reader, "chr1")
coords<-getRefCoords(reader, "chr1")
rg<-bamRange(reader,coords)
bamClose(reader)
}
\keyword{getRefData}
