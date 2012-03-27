\name{cigarSize}
\alias{cigarSize}
\title{Returns size of cigarData in bamAlign}
\description{
Size of cigarData denotes the number of cigar-items (e.g. 3 for cigarData 45M3501N56M).
}
\usage{cigarSize(object)}

\arguments{
  \item{object}{ bamAlign object }
}
\value{Integer.}
\author{Wolfgang Kaisers}

\examples{
library(rbamtools)
bam_file <- system.file("extdata", "accepted_hits.bam", package="rbamtools")
reader<-bamReader(bam_file)
align<-getNextAlign(reader)
name(align)
position(align)
cigarSize(align)
bamClose(reader)
}
\keyword{cigarSize}
\keyword{cigarData}
\keyword{bamAlign}
