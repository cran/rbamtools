\name{GenomePartition-class}
\Rdversion{1.1}
\docType{class}
%
\alias{GenomePartition-class}
%
\alias{countPartition}
\alias{countPartition-methods}
\alias{countPartition,GenomePartition,bamReader-method}
\alias{countPartition,GenomePartition,data.frame-method}
%
\alias{checkPartition}
\alias{checkPartition-methods}
\alias{checkPartition,GenomePartition,data.frame-method}
%
\alias{genomePartition}
\alias{genomePartition-methods}
\alias{genomePartition,bamReader,data.frame-method}
%
\alias{getRefData,GenomePartition-method}
%
\alias{getFileTable}
\alias{getFileTable-methods}
\alias{getFileTable,GenomePartition-method}
%
\alias{getSeqNr}
\alias{getSeqNr-methods}
\alias{getSeqNr,GenomePartition-method}
%
\alias{getAlignCounts}
\alias{getAlignCounts-methods}
\alias{getAlignCounts,GenomePartition-method}
%
\alias{getGridAlignCounts}
\alias{getGridAlignCounts-methods}
\alias{getGridAlignCounts,GenomePartition-method}
%
\alias{show,GenomePartition-method}
%
\title{Class \code{"GenomePartition"}: .}
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
% Description
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
\description{
.
}
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
% Creation of objects
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
\section{Objects from the Class}{can be created using the 
\code{genomePartition} function on an opened bamReader and a data.frame
containing target regions (usually gene annotation data).
}
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
% Slots
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
\section{Slots}{
    \describe{
        \item{\code{ev}:}{\code{"environment"}. Environment which
            contains data as data.frames.}
  }
}
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
% Details
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
\details{The environment contains at least three objects:
\code{reftable}, \code{reflist} and \code{genome}. The reftable data.frame
contains data on reference sequences (chromosomes) derived from the
bamReader \code{getRefData} function. The reflist is a list containing
calculated partition segments based on the \code{genome} argument of the 
genomePartition function and a underlying equidistant grid on the genome.
After counting BAM alignments using a genomePartition object, the
environment additionally contains a filetable data.frame.
}
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
% Methods
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
\section{Methods}{
    \describe{
        \item{show}{\code{signature(object="GenomePartition")}:
            Prints a short message with some summarizing data.}
        %
        \item{getSeqNr}{\code{signature(object="GenomePartition")}:
            Returns number of reference sequences in refdata.}
        %
        \item{getRefData}{\code{signature(object="GenomePartition")}:
            Returns refdata data.frame.}
        %
        \item{getPos}{\code{signature(object="GenomePartition")}:
            Returns numeric position values for align depth.}
        %
        \item{countPartition}{\code{signature(partition="GenomePartition", src="bamReader")}:
            Counts BAM alignments on the contained grid from a single
            BAM file (represented by an opened bamReader.}
        %
        \item{countPartition}{\code{signature(partition="GenomePartition", src="data.frame")}:
            Counts BAM alignments on the contained grid from multiple
            BAM files (filenames and sample names given in src)}
        %
        \item{checkPartition}{\code{signature(partition="GenomePartition", src="data.frame", verbose="logical")}:
            Checks incoming arguments and present reference sequences
            for consistency. Intended to be used as preparation for
            \code{countPartition} in order to prevent unexpected 
            routine terminations.}
        %
        \item{getFileTable}{\code{signature(object="GenomePartition")}:
            Returns the contained filetable (constituted by usage of
            \code{countPartition} function).}
        %
        \item{getAlignCounts}{\code{signature(object="GenomePartition")}:
            Returns counted alignments in annotated regions (exons, genes).}
        %
        \item{getGridAlignCounts}{\code{signature(object="GenomePartition")}:
            Returns counted alignments in underlying equidistant grid.}
}}
\author{Wolfgang Kaisers}
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
% Examples
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
\examples{
# Open (indexed) BAM file
bam<-system.file("extdata", "accepted_hits.bam", package="rbamtools")
reader<-bamReader(bam,idx=TRUE)

# Provide exon positions
id <- 1:13
seqid <- "chr1"
gene <- "WASH7P"
ensg_id <- "ENSG00000227232"
start <- c(14411, 15000, 15796, 15904, 16607, 16748, 16858, 17233,
           17602, 17915, 18268, 24737, 29534)
end <-   c(14502, 15038, 15901, 15947, 16745, 16765, 17055, 17364,
           17742, 18061, 18366, 24891, 29806)
           
ref <- data.frame(id=id, seqid=seqid, begin=start, end=end, gene=gene, ensg=ensg_id)

# Create partition (adds equidistant grid)
partition <- genomePartition(reader, ref)

# Count alignments
countPartition(partition, reader)

# Extract result data
urc <- getAlignCounts(partition)
gac <- getGridAlignCounts(partition)

}
\keyword{bamAlign}
