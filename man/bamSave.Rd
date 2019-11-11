\name{bamSave}
%
\alias{bamSave}
\alias{bamSave-methods}
%\alias{bamSave,bamReader-method}
\alias{bamSave,bamReader,ANY,ANY-method}
\alias{bamSave,bamReader,bamWriter,missing-method}
\alias{bamSave,bamWriter,bamAlign,numeric-method}
\alias{bamSave,bamWriter,bamRange,character-method}
\alias{bamSave,bamWriter,bamRange,missing-method}
\alias{bamSave,bamWriter,bamRange,numeric-method}
%
\title{bamSave: Creation of a sorted copy of a BAM file.}
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
% Description
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
\description{The methods for \code{bamSave} write BAM alignment data
(either single alignments as \code{bamAlign} objects or 
lists of alignments as \code{bamRange} objects) out to a bamWriter.

}
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
% Details
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
\details{
The different methods specialize on different methods for determining
the Reference sequence ID in the output file.
}
%"bamWriter", "bamRange", "character"
\usage{bamSave(object, value, refid, ...)}
\arguments{
    \item{object}{\code{bamWriter}. The writer must be opened.}
    %
    \item{value}{Object of type \code{bamAlign} or \code{bamRange}.}
    %
    \item{refid}{(Optional) Defines sequence association (ID) in 
        Reference sequence dictionary of output BAM file.
        May be given as numeric, or character. When missing, ID
        is determined from reference name (only for \code{bamRange}).}
    %
    \item{...}{Currently unused.}

}
\author{Wolfgang Kaisers}
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
% Examples
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
\examples{
bam<-system.file("extdata","accepted_hits.bam",package="rbamtools")
reader<-bamReader(bam, idx=TRUE)
rg <- bamRange(reader, getRefCoords(reader,"chr1"))
\donttest{
writer <- bamWriter(bamReader, "out.BAM")
bamSave(writer, rg)
bamClose(writer)
}

}
\keyword{bamSave}
\keyword{bamReader}
