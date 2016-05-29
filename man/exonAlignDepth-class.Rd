\name{exonAlignDepth-class}
\Rdversion{1.1}
\docType{class}
%
\alias{exonAlignDepth-class}
\alias{initialize, exonAlignDepth}
%
\alias{[,exonAlignDepth,ANY,ANY,ANY-method}
%
\alias{show,exonAlignDepth-method}
\alias{length,exonAlignDepth-method}
%
\alias{exonAlignDepth}
\alias{exonAlignDepth-methods}
\alias{exonAlignDepth,geneAlignDepth-method}
%
\alias{aldRatio}
\alias{aldRatio-methods}
\alias{aldRatio,exonAlignDepth-method}
%
\alias{cutFlatAlignDepth}
\alias{cutFlatAlignDepth-methods}
\alias{cutFlatAlignDepth,exonAlignDepth-method}
%
\alias{getNormFactor}
\alias{getNormFactor-methods}
\alias{getNormFactor,exonAlignDepth-method}
%
\alias{groupAldMatrix}
\alias{groupAldMatrix-methods}
\alias{groupAldMatrix,exonAlignDepth-method}
%
\alias{groupAldTable}
\alias{groupAldTable-methods}
\alias{groupAldTable,exonAlignDepth-method}
%
\alias{junctionSites}
\alias{junctionSites-methods}
\alias{junctionSites,exonAlignDepth-method}
%
\title{Class \code{"exonAlignDepth"}: Align-depth data for single genes
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
Objects can be created by calls of the form \code{exonAlignDepth(sal)}.
}
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
% Slots
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
\section{Slots}{
    \describe{
        \item{ald}{\code{"matrix"}:
            Matrix containing alignment depth data.}
        %
        \item{aldRatio}{\code{"data.frame"}:
            Contains mean alignment depth and
            alignment-depth-ratio data.
            \code{nr} is the ratio between the in-place and
            next alignment depth.
            \code{pr} is the ratio between
            in-place and previous alignment depth.}
        %
        \item{junctions}{\code{"data.frame"}:
            Contains alignment-gap-sites which had been generated
            when reading \code{geneAlignDepth} data.}
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
        \item{\code{nAligns}:}{\code{"numeric"}:
            Total number of alignments in each BAM file.}
        %
        \item{\code{group}:}{\code{"factor"}:
            Group assignment}
        %
        \item{label}{\code{"character"}:
            Short textual identifier for each sample.}
        %
}
}
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
% Methods
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
\section{Methods}{
    \describe{
        \item{show}{\code{signature(object="exonAlignDepth")}:
            Prints a short message with some summarizing data.}
        %
        \item{exonAlignDepth}{\code{signature(object=c("sampleBamFiles", "geneModel"}:
            Constructs \code{exonAlignDepth} object and
            reads align depth data from BAM files.}
        %
        \item{plot}{\code{signature(object="exonAlignDepth")}:
            Plots align depth data.}
        %
        \item{aldRatio}{\code{signature(object="exonAlignDepth")}:
            Returns data.frame containing alignment-depth-ratio data
            used for identification of exon junctions.}
        %
        \item{junctionSites}{\code{signature(object="exonAlignDepth")}:
            Returns data.frame containing junction positions.}
        %
        \item{getNormFactor}{\code{signature(object="exonAlignDepth")}:
            Returns numeric value which will be uses for
            plotting and \code{groupAldMatrix} and \code{groupAldTable}.}
        %
        \item{groupAldMatrix}{\code{signature(object="exonAlignDepth")}:
            Returns matrix containing mean alignment depth values.
            Data for each sample is stored in one column.
            Data for each genomic position is stored in one row. A
            summarizing function \code{f} may be given
            (Default is \code{mean}).}
        %
        \item{groupAldTable}{\code{signature(object="exonAlignDepth")}:
            Returns data.frame containing three columns.
            The first column contains genomic positions,
            the second position contains group assignment and
            the third position contains (normalized) alignment
            dept values. A
            summarizing function \code{f} may be given
            (Default is \code{mean}).}
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
## - - - - - - - - - - - - - - - - - - - - - - ##
# Extract exonAlignDepth object
ead <- exonAlignDepth(gad, ratioLim=5, infVal=1000)
ald <- aldRatio(ead)
jc <- junctionSites(ead)
getNormFactor(ead)
cead <- cutFlatAlignDepth(ead, ratio=0.1)
}
\keyword{bamAlign}
