\name{countNucs}
\Rdversion{1.1}
\docType{class}
\alias{countNucs}
\alias{countNucs-methods}
\alias{countNucs,bamAlign-method}
\alias{countNucs,bamRange-method}
\title{Counting nucleotides in \code{bamAlign} and \code{bamRange}}
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
% Description
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
\description{The function counts occurrence of the Nucleotides A,C,G,T in 
\code{bamAlign} and \code{bamRange} objects. 
Any other values will be combined counted in the last value (N).
The function returns an integer vector of length 5. 
The names indicate which position contains the count value each nucleotide.
}
\usage{countNucs(object)}
\arguments{
  \item{object}{\code{bamAlign} or \code{bamRange} object.}
}
\value{Integer (of length 4).}
\author{Wolfgang Kaisers}
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
% Examples
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
\examples{
# A) For bamAlign
align<-bamAlign("HWUSI-0001","ACCGGGTTTT","Qual/Strng","4M10N6M",refid=0,position=100)
countNucs(align)

# B) For bamRange
bam<-system.file("extdata","accepted_hits.bam",package="rbamtools")
reader<-bamReader(bam,idx=TRUE)
coords<-c(0,0,14730)
range<-bamRange(reader,coords)
countNucs(range)
}
\keyword{countNucs}
\keyword{bamAlign}
\keyword{bamRange}
