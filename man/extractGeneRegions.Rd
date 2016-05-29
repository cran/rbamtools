\name{extractGeneRegions}
%
\alias{extractGeneRegions}
\alias{extractGeneRegions-methods}
\alias{extractGeneRegions,bamReader,bamWriter,geneList-method}
\alias{extractGeneRegions,bamReader,character,geneList-method}
\alias{extractGeneRegions,sampleBamFiles,sampleBamFiles,geneList-method}
%
\title{Extraction of alignments from given genetic regions
and BAM files into a set of output BAM files.}
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
% Description
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
\description{The function \code{extractGeneRegions} writes aligns or
some part of the stored aligns in a BAM file as reads into a
BAM output file.}
\usage{extractGeneRegions(src, trg, gl)}
\arguments{
    \item{src}{\code{sampleBamFiles}: List of source BAM files where alignments are
        read from. BAM index files must exist and be present in src object.}
    %
    \item{trg}{\code{sampleBamFiles}: List of target BAM files where alignments are
        written to.}
    %
    \item{gl}{\code{geneList}: Objects from which genetic regions are taken.}
}
\value{Numeric. Number of written alignments for each file.}
%
\details{Reference sequence names in \code{geneList} object must all match
sequence names in source BAM files.}
%
\author{Wolfgang Kaisers}
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
% Examples
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
\examples{
# + + + + + + + + + + + + + + + + + + #
# A) Create sampleBamFiles objects
# + + + + + + + + + + + + + + + + + + #
bam <- system.file("extdata","accepted_hits.bam",package="rbamtools")
outdir <- tempdir()
# Input sampleBamFiles object
ibs <- sampleBamFiles(bam)
# Output sampleBamFiles object
obs <- sampleBamFiles(file.path(outdir, "outu.bam"))
# + + + + + + + + + + + + + + + + + + #
# B) Create geneList Object
# + + + + + + + + + + + + + + + + + + #
rfile <- system.file("extdata",
            "hs.ucsc.small.RData",
                package="refGenome")

ucr <- loadGenome(rfile)
gt <- getGeneTable(ucr)
gl <- geneList(ucr, gt$gene_id)
\dontrun{
extractGeneRegions(ibs, obs, gl)
bamFiles(obs) <- file.path(outdir, bamSort(obs, "out"))
bamIdxFiles(obs) <- paste(bamFiles(obs), "bai", sep=".")
createIndex(obs)
}
}
\keyword{extractGeneRegions}
\keyword{extractGeneRegions}
