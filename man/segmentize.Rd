\name{segmentize}
%
\Rdversion{1.1}
%
\alias{segmentize}
\alias{segmentize-methods}
\alias{segmentize,ANY-method}
\alias{segmentize,matrix-method}
\alias{segmentize,data.frame-method}
%
\title{Segmentation of vector using segment coordinates}
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
% Description
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
\description{The function returns extracts of a data vector using
two vector defining begin and end coordinates. An offset can be used.
}
\usage{segmentize(x, begin, end, offset=1, margin=1, invert=FALSE)}
\arguments{
    \item{x}{\code{numeric} object.
        Data vector from which segments are extracted.}
    %
    \item{begin}{\code{numeric} object.
        Indices first segment elements.}
    %
    \item{end}{\code{numeric} object.
        Indices of last segment elements. Length must be equal
        to length(begin). \code{segmentize} requires all(end>=begin).}
    %
    \item{offset}{\code{numeric} object.
        The indices are shifted so that the first entry is indexed
        by value indicated by \code{offset}. 
        Standard value is 1 (default behaviour in R).}
    %
    \item{margin}{\code{numeric} object.
        Direction of segmentation for \code{matrix} and \code{data.frame}
        objects (1=rows, 2=columns).}
    %
    \item{invert}{\code{logical} object.
        When \code{TRUE}, the complement of the selected segments is returned.
        Here, the ordering of the returned values cannot be changed.}
}
\details{The inversion option internally operates using logical values.
As a consequence, when \code{begin} and \code{end} define overlapping
regions, in effect the union of these regions is removed,
because regions are only retained when they are not covered by any
segment defined by \code{begin} and \code{end}.
}
\value{Integer vector}
\author{Wolfgang Kaisers}
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
% Examples
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
\examples{
# Create data vector
x <- rep(0, 11)
x[3:5] <- 1:3
x[7:9] <- 4:6
names(x) <- 10:20
# Define extracted segments
sgb <- c(3, 7)
sge <- c(5, 9)
sgm <- segmentize(x, sgb, sge)
names(sgm) <- segmentize(names(x), sgb, sge)
# Use offset
offset <- 10
sgb <- c(12, 16)
sge <- c(14, 18)
segmentize(x, sgb, sge, offset)
# Matrix
m <- matrix(0L, nrow=11, ncol=5)
rownames(m) <- 10:20
colnames(m) <- letters[1:5]

for(i in 1:5)
    m[c(3:5, 7:9), i] <- c(1:3, 4:6) * i
ms <- segmentize(m, sgb, sge, offset)
}
\keyword{segmentize}
