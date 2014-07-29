\name{nucStats}
\Rdversion{1.1}
\docType{class}
\alias{nucStats}
\alias{nucStats-methods}
\alias{nucStats,bamReader-method}
\alias{nucStats,character-method}
\alias{nucStats,bamRange-method}
\title{Table nucleotides in whole BAM file(s)}
\description{The function counts occurrence of the Nucleotides A,C,G,T for whole BAM file from \code{bamReader} object or a list of BAM files. Letters different from A,C,G,T are subsumed in the value N. The function returns a \code{data.frame} with 8 columns. The names indicate which position contains the count value each nucleotide. The function also calculates GC-content and AT/GC ratio which make up the last two columns. The returned data.frame contains one row for each reference sequence (when used for \code{bamReader}) or one row for each given BAM file name.
}
\usage{nucStats(object,...)}
\arguments{
  \item{object}{\code{bamReader} object. The reader must be opened and have initialized index.}
  \item{...}{Optional name for BAM-index files (argument name=idxInfiles).}
}
\value{data.frame}
\author{Wolfgang Kaisers}
\examples{
bam<-system.file("extdata","accepted_hits.bam",package="rbamtools")
nucStats(bam)
idx<- system.file("extdata", "accepted_hits.bam.bai", package="rbamtools")
nucStats(bam,idxInfiles=idx)
reader<-bamReader(bam,idx=TRUE)
nucStats(reader)
range<-bamRange(reader,c(0,0,14000))
nucStats(range)
}
\keyword{nucStats}
\keyword{bamReader}
\keyword{GC content}
\keyword{AT/CG ratio}
