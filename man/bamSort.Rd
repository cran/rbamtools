\name{bamSort}
\alias{bamSort}
\title{bamSort(bamAlign)}
\description{Sorts a given BAM-file (by position or by name) into a second BAM-file.}
\usage{bamSort(object,prefix,byName=FALSE,maxmem=100000000)}
\arguments{
  \item{object}{An Instance of class bamReader}
  \item{prefix}{Character String}
  \item{byName}{Logical}
  \item{maxmem}{Approximate maximum memory usage}
}
\details{\verb{
BAM-files must be sorted before an index can be created.
During this routine, some content of the source file
is copied into the working memory, sorted and written
into temporary files. Finally, the sorted files are
merged into a single target. The name of the target
is [prefix].bam.
Number and size of temporary files depend on source
file size and used working memory (maxmem). Small
maxmem values result in large numbers of temporary
files.
The minimum maxmem value is 100000000. Smaller maxmem
values cause an error.
}}
\value{None}
\author{Wolfgang Kaisers}
\examples{
library(rbamtools)
bam_file <- system.file("extdata", "accepted_hits.bam", package="rbamtools")
reader<-bamReader(bam_file)
isOpen(reader)
bamSort(reader,"sorted")
bamClose(reader)
}
\keyword{bamReader}
\keyword{bamSort}
