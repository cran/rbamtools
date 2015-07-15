\name{extractBamRegions}
\alias{extractBamRegions}
\title{extractBamRegions: Extraction of alignments from given (genetic) regions
and BAM files into a set of output BAM files..}
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
% Description
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
\description{The function \code{extractBamRegions} writes (remaining) aligns or
some part of the stored aligns in a BAM file as reads into a
compressed fastq file.}
\usage{extractBamRegions(bamFiles, ranges, idxFiles, outFiles)}
\arguments{
    \item{bamFiles}{\code{character}: Vector of BAM file names. }
    %
    \item{ranges}{\code{data.frame}: Must contain three columns 
        \code{seqid}, \code{start} and \code{end}. The seqid }
    %
    \item{idxFiles}{\code{character}: Vector of BAM index file names.}
    %
    \item{outFiles}{\code{character}: Vector of output file names.}
}
\value{No return value given}
%
\details{bamFiles, idxFiles and outFiles are checked for equal length.
bamFiles and idxFiles are checked for existance. 
There are sensible default values given for idxFiles and outFiles. The
default location for writing output files is the current working
directory.}
%
\author{Wolfgang Kaisers}
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
% Examples
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
\examples{
bam<-system.file("extdata","accepted_hits.bam",package="rbamtools")
gene_pos <- data.frame(seqid="chr1", start=15e3, end=20e3)
\dontrun{extractBamRegions(bam, gene_pos)}
}
\keyword{extractBamRegions}
\keyword{extractBamRegions}
