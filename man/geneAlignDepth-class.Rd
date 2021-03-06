\name{geneAlignDepth-class}
\Rdversion{1.1}
\docType{class}
%
\alias{geneAlignDepth-class}
\alias{initialize, geneAlignDepth}
%
\alias{[,geneAlignDepth,ANY,ANY,ANY-method}
%
\alias{show,geneAlignDepth-method}
\alias{length,geneAlignDepth-method}
%
\alias{geneAlignDepth}
\alias{geneAlignDepth-methods}
\alias{geneAlignDepth,sampleBamFiles,geneModel-method}
%
\title{Class \code{"geneAlignDepth"}: Align-depth data for single genes
on multiple samples (BAM files).}
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
% Description
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
\description{Objects of this class combine data from \code{geneModel}
and from \code{sampleBamFiles} objects allowing to read BAM alings for
gene defined regions from multiple BAM files.
\code{geneAlignDepth} object contain an Align depth matrix
which can be visualized using a generic plot function.
}
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
% Creation of objects
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
\section{Objects from the Class}{
Objects can be created by calls of the form \code{getgeneAlignDepth(gesa, gm)
}.
}
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
% Slots
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
\section{Slots}{
    \describe{
        \item{\code{bamFiles}:}{\code{"character"}:
            Location of BAM files}
        %
        \item{\code{bamIdxFiles}:}{\code{"character"}:
            Location of BAM index files}
        \item{\code{nAligns}:}{\code{"numeric"}:
            Total number of alignments in each BAM file.}
        %
        \item{\code{group}:}{\code{"factor"}:
            Group assignment}
        %
        \item{label}{\code{"character"}:
            Short textual identifier for each sample.}
        %
        \item{length}{\code{"integer"}:
            Vector length for bamFiles, bamIdxFiles, nAligns,
            group and label.}
        %
        \item{ald}{\code{"matrix"}:
            Matrix containing alignment depth data.}
        %
        \item{gapSites}{\code{"data.frame"}:
            Align gap sites data from genetic region.}
        %
        \item{ev}{\code{"environment"}:
            Contains additional data  (e.g. group table).}
        %
        \item{\code{gene_id}:}{\code{"character"}:
            Gene identifier (e.g. Ensembl or UCSC).}
        %
        \item{\code{gene_name}:}{\code{"character"}:
            Gene name.}
        %
        \item{\code{seq_name}:}{\code{"character"}:
            Reference sequence (chromosome) name.}
        %
        \item{\code{strand}:}{\code{"character"}:
            Strand orientation of gene on reference sequence (+, - or *)}
        %
        \item{\code{coords}:}{\code{"numeric"}:
            Gene coordinates (i.e. start and end position).}
        %
  }
}
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
% Methods
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
\section{Methods}{
    \describe{
        \item{show}{\code{signature(object="geneAlignDepth")}:
            Prints a short message with some summarizing data.}
        %
        \item{geneAlignDepth}{\code{signature(object=c("sampleBamFiles", "geneModel"}:
            Constructs \code{geneAlignDepth} object and
            reads align depth data from BAM files.}
        %
        \item{plot}{\code{signature(object="geneAlignDepth")}:
            Plots align depth data.}
        %
}}
\author{Wolfgang Kaisers}
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
% Examples
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
\examples{
## - - - - - - - - - - - - - - - - - - - - - - ##
# Construct sampleBamFiles object
bam<-system.file("extdata", "accepted_hits.bam", package="rbamtools")
bs <-sampleBamFiles(1)
bamFiles(bs) <- bam
sampleLabels(bs) <- "s1"
sampleGroups(bs) <- "g1"
checkBamFiles(bs)
nAligns(bs) <- bamCountAll(bs)
bs
## - - - - - - - - - - - - - - - - - - - - - - ##
# Construct geneModel object
library(refGenome)
ucfile<-system.file("extdata", "hs.ucsc.small.RData", package="refGenome")
uc<-loadGenome(ucfile)
gt <- getGeneTable(uc)
gene_id <- as.character(gt$gene_id[1])
gm <- geneModel(uc, gene_id)
## - - - - - - - - - - - - - - - - - - - - - - ##
# Construct geneAlignDepth object
gad <- geneAlignDepth(bs, gm)
plot(gad, col="gray50")
}
\keyword{bamAlign}
