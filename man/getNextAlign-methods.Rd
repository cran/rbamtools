\name{getNextAlign-methods}
\docType{methods}
\alias{getNextAlign}
\alias{getNextAlign-methods}
\title{\code{getNextAlign}: Retrieving next align from bamReader or bamRange}
\description{
Returns an object of class bamAlign from bamReader or bamRange.
}
\section{Methods}{
\describe{
\item{\code{signature(object = "bamRange")}}{
An instance of class bamRange.
}
\item{\code{signature(object = "bamReader")}}{
An instance of class bamReader
}
}}
\value{An object of class bamAlign or NULL.}
\author{Wolfgang Kaisers}
\examples{
library(rbamtools)
bam<-system.file("extdata","accepted_hits.bam",package="rbamtools")
idx<-system.file("extdata","accepted_hits.bam.bai",package="rbamtools")
reader<-bamReader(bam)
isOpen(reader)
align<-getNextAlign(reader)
if(is.null(align))
	print("End of File!\n")
load.index(reader,idx)
index.initialized(reader)
coords<-as.integer(c(0,0,249250621))
range<-bamRange(reader,coords)
align<-getNextAlign(range)
position(align)
bamClose(reader)
}
\keyword{methods}
\keyword{getNextAlign}
\keyword{bamReader}
\keyword{bamRange}
