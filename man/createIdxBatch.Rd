\name{createIdxBatch}
\alias{createIdxBatch}
\alias{create.idx.batch}
\title{createIdxBatch: Creation of index files for multiple BAM files.}
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
% Description
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
\description{The function takes a vector of BAM-file names 
(plus optionally accompanying names of BAM-index-files) and checks for 
existing BAM-index files.
When index files do not exist, the function creates the missing index files.}
\usage{createIdxBatch(bam,idx=paste(bam,".bai",sep=""),rebuild=FALSE)}
\arguments{
    \item{bam}{Filenames of BAM-files to be opened.}
    \item{idx}{Optional: Name of BAM-index files}
    \item{rebuild}{Optional: Logical value. 
        When \code{TRUE} the function rebuilds existing BAM-index files.}
}
\author{Wolfgang Kaisers}
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
% Examples
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
\examples{
\dontrun{bam<-system.file("extdata","accepted_hits.bam",package="rbamtools")}
\dontrun{createIdxBatch(bam)}
}
\keyword{createIdxBatch}
