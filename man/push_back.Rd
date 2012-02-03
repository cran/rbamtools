\name{push_back}
\alias{push_back}
\alias{push_front}
\alias{pop_back}
\alias{pop_front}
\title{Pushing and popping bamAligns on or off a bamRange list.}
\description{A bamRange is a (linked-) list of bamAlign objects. The four functions add and remove single bamAlings at the beginning or the end of the list.}
\usage{
push_back(object, value)
push_front(object, value)
pop_back(object)
pop_front(object)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{A bamRange object.}
  \item{value}{A bamAlign object.}
}
\details{
push_back adds the given bamAlign as last item in the list. pop_back removes the last item from the bamRange list. push_fron adds the given bamAlign as first item in the list. pop_front removes the first item from the bamRange list. The pop-functions have no effect when the bamRang list is empty.
}
\value{None.}
\author{Wolfgang Kaisers}
\examples{
bam_file <- system.file("extdata", "accepted_hits.bam", package="rbamtools")
idx_file <- system.file("extdata", "accepted_hits.bam.bai", package="rbamtools")
reader<-bamReader(bam_file)
isOpen(reader)
loadIndex(reader,idx_file)
index.initialized(reader)
coords<-as.integer(c(0,0,1874))
range<-bamRange(reader,coords)
align<-getNextAlign(range)
push_back(range,align)
pop_back(range)
push_front(range,align)
pop_front(range)
bamClose(reader)
}
\keyword{push_back}
\keyword{push_front}
\keyword{pop_back}
\keyword{pop_front}
