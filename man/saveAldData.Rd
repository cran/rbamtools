\name{saveAldData}
\alias{saveAldData}
\alias{saveAldData-methods}
\alias{saveAldData,sampleBamFiles,geneList-method}
\title{saveAldData: Save alignment depth plots for single genes}
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
% Description
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %}
\description{The function \code{saveAldData} reads alignment data for
    all genes in \code{geneList} sequentially from all Bam files
    listed in \code{sampleBamFiles}.
    The function cuts out intronic regions using the \code{exonAlignDepth}
    class and calculates a smoothed groupwise estimate
    using the \code{exonLoessModel} function.
    In a second region cutting step, remaining regions with low alignment
    depth are cut out using \code{cutFlatAlignDepth}.
    Finally, the function calculates group-wise ratios of alignment
    depth using the \code{groupRatio} function.

    }
\usage{saveAldData(bs, gl, path, order=NULL, lim=1.1, startId=1, f=mean)}
\arguments{
    \item{bs}{sampleBamFiles. Contains BAM file names and group assignment.}
    %
    \item{gl}{geneList. Contains a list of \code{geneModel} objects.
        Only gene-name, gene-id and global genetic region is used.
        Object construction can be abridged using \code{interior=FALSE}.}
    %
    \item{path}{character. Base directory where output files are written
    to.}
    %
    \item{order}{numeric. Defines order in which group-levels of
        \code{sampleBamFiles} object are used for \code{groupRatio}}
    %
    \item{lim}{numeric. Limit used for \code{groupRatio} function.
        The percentage of genomic position where this limit is
        exceeded between all (ordered) group levels is
        given as gr (and cgr) value in \code{aldrat} table.}
    %
    \item{startId}{numeric. Index used for first gene. Starting
        value for subsequent numbering of genes.}
    %
    \item{f}{function. Function for calculation of group accumulates.
        Defaults to \code{mean}. Alternatively \code{median} may also
        be used.
    }
}
\author{Wolfgang Kaisers}
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
% Examples
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
\examples{
bam <- system.file("extdata", "accepted_hits.bam", package="rbamtools")
reader <- bamReader(bam, idx=TRUE)
# Extract data for HLHL17 gene:
seqid <- "chr1"
start <- 895967
end <- 901099
gene_name <- "HLHL17" # (optional)
ranges <- data.frame(seqid=seqid, start=start, end=end, gene_name=gene_name)
\dontrun{saveAldData(reader, ranges=ranges, filename="new_file.bam", f=mean)}
}
\keyword{saveAldData}
\seealso{bamReader}
