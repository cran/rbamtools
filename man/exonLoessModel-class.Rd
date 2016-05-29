\name{exonLoessModel-class}
\Rdversion{1.1}
\docType{class}
%
\alias{exonLoessModel-class}
\alias{initialize, exonLoessModel}
%
\alias{[,exonLoessModel,ANY,ANY,ANY-method}
%
\alias{show,exonLoessModel-method}
\alias{length,exonLoessModel-method}
%
\alias{exonLoessModel}
\alias{exonLoessModel-methods}
\alias{exonLoessModel,exonAlignDepth-method}
%
\alias{cutFlatAlignDepth,exonLoessModel-method}
%
\title{Class \code{"exonLoessModel"}: Align-depth data for single genes
on multiple samples (BAM files). Intronic regions are cut out
of align depth matrix.}
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
% Description
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
\description{Objects of this class combine data from \code{geneAlignDepth}.
During construction, the align depth matrix is segmentized, so that
only (unified) exonic regions are represented in align depth matrix.
}
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
% Creation of objects
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
\section{Objects from the Class}{
Objects can be created by calls of the form \code{exonLoessModel(sal)}.
}
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
% Slots
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
\section{Slots}{
    \describe{
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
        \item{\code{nAligns}:}{\code{"numeric"}:
            Total number of alignments in each BAM file.}
        %
        \item{\code{group}:}{\code{"character"}:
            Group assignment}
        %
        \item{\code{loessPred}:}{\code{"matrix"}:
            Loess estimated align depth values for each group.
            Each column contains estimates for one group.}
}
}
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
% Methods
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
\section{Methods}{
    \describe{
        \item{show}{\code{signature(object="exonLoessModel")}:
            Prints a short message with some summarizing data.}
        %
        \item{exonLoessModel}{\code{signature(object=c("sampleBamFiles", "geneModel"}:
            Constructs \code{exonLoessModel} object and
            reads align depth data from BAM files.}
        %
        \item{plot}{\code{signature(object="exonLoessModel")}:
            Plots align depth data.}
        %
        \item{groupRatio}{\code{signature(object="exonLoessModel", lim="numeric")}:
            Returns single numeric value. From align depth estimates in
            loessPred matrix, the ratio between a column and the
            precedent column is calculated. For each position, the
            minimum ratio is extracted. The function returns the
            relative number of positions, where the minimum ratio
            exceeds given \code{lim} value.}
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
## - - - - - - - - - - - - - - - - - - - - - - ##
# Extract exonLoessModel object
ead <- exonAlignDepth(gad, ratioLim=5, infVal=1000)
elm <- exonLoessModel(ead)
celm <- cutFlatAlignDepth(elm, ratio=0.1)
groupRatio(celm, lim=1.2, cut=0, order=1)
}
\keyword{bamAlign}
