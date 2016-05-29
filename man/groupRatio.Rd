\name{groupRatio}
%
\alias{groupRatio}
\alias{groupRatio-methods}
\alias{groupRatio,exonLoessModel-method}
%
\title{Calculates group-wise ratios of alignment depth (AD)}
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
% Description
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
\description{\code{groupRatio} takes a \code{bamRange} and returns a data.frame
(128 rows, number of columns=length of the longest sequence in range).
\code{groupRatio} counts occurrences for every sequence position (column) and
every phred value (row).
\code{getQualQuantiles} takes a \code{bamReader} and a vector of quantiles
(must be between 0 and 1) and returns a data.frame.
The data.frame contains one row for each quantile and also as many columns
as the maximum sequence length. \code{plotQualQuant} plots the values
for quanties 0.1,0.25,0.5,0.75 and 0.9.}
%
\usage{groupRatio(object, lim=1.2, cut=0, order=NULL, f=mean)}
\arguments{
    \item{object}{exonLoessModel}
    %
    \item{lim}{numeric. Limit ratio. Must be > 1. The function returns
    the fraction of genetic position where AD-ratio between groups is > lim
    or the fraction of positions where AD-Ratio is < 1/lim
    (i.e the larger ratio).}
    %
    \item{cut}{numeric. When >0 , the function uses \code{cutFlatAlignDepth}
    for cutting out low alignment depth regions before calculating.
    alignment depth ratio.}
    %
    \item{order}{numeric. When given, the function reorders the sample groups.
    Can be used to provide ascending (or descending) group ordering, e.g.
    group1 < group2 < group3.}
    %
    \item{f}{function. Function for calculation of group accumulates.
        Defaults to \code{mean}. Alternatively \code{median} may also
        be used.
    }
    %
}
\value{numeric}
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
% Details
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
\details{The size of the returned value (abs(groupRatio)) indicates on which
proportion of the genetic region, AD ratio between subsequent
groups exceeds the given limit. For lim=1.1, group1<group2<group3, a
returned value of 0.8 says that the AD ratios group2:group1 and group3:group2
are at least 1.1 (> 1) on 80 percent of the contained genomic positions.
Negative values say that the relation is group1>group2>group3. This
allows discrimination of up- and down-regulated genes.}
\author{Wolfgang Kaisers}
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
% Examples
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
\examples{
## - - - - - - - - - - - - - - - - - - - - - - ##
# Construct sampleBamFiles object
bam <- system.file("extdata", "accepted_hits.bam", package="rbamtools")
bs <- sampleBamFiles(1)
bamFiles(bs) <- bam
sampleLabels(bs) <- "s1"
sampleGroups(bs) <- "g1"
checkBamFiles(bs)
nAligns(bs) <- bamCountAll(bs)
bs
## - - - - - - - - - - - - - - - - - - - - - - ##
# Construct geneModel object
library(refGenome)
ucfile <- system.file("extdata", "hs.ucsc.small.RData", package="refGenome")
uc <- loadGenome(ucfile)
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
\keyword{groupRatio}
\keyword{bamRange}
