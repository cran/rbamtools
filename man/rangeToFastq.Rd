\name{rangeToFastq}
%
\alias{rangeToFastq}
\alias{rangeToFastq-methods}
\alias{rangeToFastq,bamRange,character-method}
%
\title{rangeToFastq: Extract read information from alignments for given
    genomic range as fastq.}
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
% Description
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
\description{The function \code{rangeToFastq} writes all (or selected) aligns
    from a \code{bamRange} into a compressed fastq file.}
\usage{rangeToFastq(object, filename, which, append=FALSE)}
\arguments{
    \item{object}{bamReader. Must be opened.}
    %
    \item{filename}{Name of output 'fastq' file.}
    %
    \item{which}{(Optional) Logical vector.
        When given, the routine checks for each of the given values
        the value of which.
        When TRUE the corresponding align will be written to the output file,
        otherwise will be skipped.}
    %
    \item{append}{(Optional) When TRUE, the routine will append to an
        existing file. Otherwise existing files will be overwritten.}
}
\author{Wolfgang Kaisers}
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
% Examples
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
\examples{
bam <- system.file("extdata", "accepted_hits.bam", package="rbamtools")
reader <- bamReader(bam, idx=TRUE)
coords <- as.integer(c(0, 0, 249250621))
range <- bamRange(reader, coords)
\donttest{rangeToFastq(range, "range.fq.gz")}
index <- sample(1:size(range), 100)
\donttest{rangeToFastq(range, "range_subset.fq.gz", which=index)}
}
\keyword{rangeToFastq}
\keyword{rangeToFastq}
