\name{countTextLines}
\alias{countTextLines}
\title{countTextLines: Counting lines in text files.}
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
% Description
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
\description{The countTextLines function takes a filename vector
and returns the number of lines contained in each file.}
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
% Details
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
\details{Fastq files usually (but not guaranteed) contain four lines of
text for each read. 
In many cases, this function can be used to count the number of
reads in a fastq file.}
\usage{countTextLines(filenames)}
\arguments{
    \item{filenames}{Character. Vector of file-names to be opened.}
}
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
% Examples
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
\examples{
filename <- system.file("extdata", "test.fastq", package="rbamtools")
countTextLines(filename)
}
\keyword{countTextLines}
\keyword{index}
