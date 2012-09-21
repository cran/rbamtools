\name{isOpen-methods}
\docType{methods}
\alias{isOpen}
\alias{isOpen-methods}
\title{ isOpen(bamReader),isOpen(bamWriter) }
\description{
Checks whether the given bamReader or bamWriter represents an opened connection to a BAM-file.
}
\usage{isOpen(con,rw="")}
\arguments{
  \item{con}{An instance of class bamReader or bamWriter.}
  \item{rw}{This argument is included solely for compatibility with the function template for connections in package base. It's not evaluated internally.}
}
\details{The function acutally checks, if the externalptr (reader or writer) slot has nil-value.}
\section{Methods}{
\describe{
\item{\code{signature(object = "bamReader")}}{ Returns TRUE when file connection is open. }
\item{\code{signature(object = "bamWriter")}}{ Returns TRUE when file connection is open. }
}}
\keyword{methods}
\keyword{isOpen}