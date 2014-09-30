\name{readPooledBamGapDf}
\alias{readPooledBamGapDf}
\title{readPooledBamGapDf: Position data for alignments gaps from multiple
    BAM files.}
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
% Description
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
\description{The function takes vectors of BAM-file names 
(plus accompanying names of BAM-index-files) and returns extracted data
as data.frame. 
The data frame contains coordinates of align gaps plus a gap-quality-score
(gqs) which quantifies information amount for detection of splice sites.}
\usage{readPooledBamGapDf(infiles, idxInfiles=paste(infiles, ".bai", sep=""))}
\arguments{
    \item{infiles}{Filenames of BAM-files to be opened for data extraction.}
    %
    \item{idxInfiles}{Optional: Name of BAM-index files}
}
\author{Wolfgang Kaisers}
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
% Examples
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
\examples{
bam <- system.file("extdata", "accepted_hits.bam", package="rbamtools")
rpb <- readPooledBamGaps(bam)
rpdf <- readPooledBamGapDf(bam)
hist(rpdf$gqs, las=1, xlab="", main="Gap quality score")
}
\keyword{readPooledBamGaps}
\keyword{readPooledBamGapDf}
