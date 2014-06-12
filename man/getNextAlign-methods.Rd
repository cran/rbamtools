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
## + + + + + + + + + + + + + + + + + + + + + + + + + + + ##
## Open bamReader
## + + + + + + + + + + + + + + + + + + + + + + + + + + + ##
bam<-system.file("extdata","accepted_hits.bam",package="rbamtools")
reader<-bamReader(bam,idx=TRUE)

## + + + + + + + + + + + + + + + + + + + + + + + + + + + ##
## Use 'getNextAlign' on bamReader
## + + + + + + + + + + + + + + + + + + + + + + + + + + + ##
align<-getNextAlign(reader)
if(is.null(align))
	print("End of File!\n")

## + + + + + + + + + + + + + + + + + + + + + + + + + + + ##
## Read bamRange
## + + + + + + + + + + + + + + + + + + + + + + + + + + + ##
coords<-as.integer(c(0,0,249250621))
range<-bamRange(reader,coords)

## + + + + + + + + + + + + + + + + + + + + + + + + + + + ##
## Use 'getNextAlign on bamRange
## + + + + + + + + + + + + + + + + + + + + + + + + + + + ##
align<-getNextAlign(range)
position(align)

bamClose(reader)
}
\keyword{methods}
\keyword{getNextAlign}
\keyword{bamReader}
\keyword{bamRange}
