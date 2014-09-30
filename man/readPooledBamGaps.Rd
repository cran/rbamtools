\name{readPooledBamGaps}
\alias{readPooledBamGaps}
\title{readPooledBamGaps: Position data for alignments gaps from multiple
    BAM files.}
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
% Description
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
\description{The function takes vectors of BAM-file names 
(plus accompanying names of BAM-index-files) and return a bamGapList object. 
The bamGapList object can be merged with other bamGapList objects and data
can be extracted with the 'as.data.frame' method.}
\usage{readPooledBamGaps(infiles,idxInfiles=paste(infiles,".bai",sep=""))}
\arguments{
    \item{infiles}{Filenames of BAM-files to be opened for data extraction.}
    %
    \item{idxInfiles}{Optional: Name of BAM-index files}
}
\value{bamGapList}
\author{Wolfgang Kaisers}
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
% Examples
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
\examples{
bam<-system.file("extdata","accepted_hits.bam",package="rbamtools")
rpb<-readPooledBamGaps(bam)
rpdf<-readPooledBamGapDf(bam)
hist(rpdf$gqs)
}
\keyword{readPooledBamGaps}
\keyword{readPooledBamGapDf}
