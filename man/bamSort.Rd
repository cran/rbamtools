\name{bamSort}
\alias{bamSort}
\title{bamSort: Creation of a sorted copy of a BAM file.}
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
% Description
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
\description{The function takes an opened instance of \code{bamReader},
a prefix for the output file and produces a sorted BAM file which is
(apart from being sorted) identical to the input file. 
BAM-files must be sorted before an index can be created.
During this routine, some content of the source file is copied into the
working memory, sorted and written into temporary files.
Finally, the sorted files are merged into a single target.
The name of the sorted BAM-file will be [prefix].bam.
Number and size of temporary files depend on source file size and used
working memory (maxmem).
Small maxmem values result in large numbers of temporary files.
The minimum maxmem value is 100000000. Smaller maxmem values cause an error.
}
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
% Details
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
\details{
The function does not take a complete name for the output file but only a
prefix.
The prefix is internally completed with a '.bam' suffix. 
This is because the samtools function 'bam_sort_core_ext' only takes a prefix.
samtools in turn produces intermediate files which also use the prefix and
which are removed again when the 'bam_sort_core_ext' cleans up.
}
\usage{bamSort(object, prefix="sorted",
            byName=FALSE, maxmem=1e+9, path=dirname(filename(object)))}
\arguments{
    \item{object}{\code{bamReader}. The reader must be opened.}
    \item{prefix}{The prefix of the output file. 
    When 'sorted' is given as prefix, the routine produces a file named
        'sorted.bam'.}
    \item{byName}{(Optional) Logical. Must have length 1.
        When TRUE the file will be sorted by align name.
        When FALSE the file will be sorted by coordinate.
        Sorting by coordinate is the prerequisite for creation of
        an index file.}
    \item{maxmem}{(Optional) Numeric. Must have length 1.
        Minimum value is 100000000. Smaller maxmem values cause an error.
        Determines how many aligns are sorted inside the working memory before
        they are written into a temporary file.}
    \item{path}{\code{character}. Directory where sorted output file
        is written to. Default value is the directory where the input
        file is located.}
}
\author{Wolfgang Kaisers}
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
% Examples
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
\examples{
bam<-system.file("extdata","accepted_hits.bam",package="rbamtools")
reader<-bamReader(bam)
bamSort(reader)
}
\keyword{bamSort}
\keyword{bamReader}
