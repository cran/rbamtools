\name{create.idx.batch}
\alias{create.idx.batch}
\title{create.idx.batch}
\description{The function takes a vector of BAM-file names (plus optionally accompanying names of BAM-index-files) and checks for existing BAM-index files. When index files do not exist, the function creates the missing index files.}
\usage{create.idx.batch(bam,idx=paste(bam,".bai",sep=""),rebuild=FALSE)}
\arguments{
  \item{bam}{Filenames of BAM-files to be opened.}
  \item{idx}{Optional: Name of BAM-index files}
  \item{rebuild}{Optional: Logical value. When \code{TRUE} the function rebuilds existing BAM-index files.}
}
\author{Wolfgang Kaisers}
\examples{
#bam<-system.file("extdata","accepted_hits.bam",package="rbamtools")
#create.idx.batch(bam)
}
\keyword{create.idx.batch}