\name{multSeq}
%
\Rdversion{1.1}
%
\alias{multSeq}
%
\title{Create combined sequences of sequential numbers.}
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
% Description
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
\description{The function takes two vectors (of type integer) of equal length.
For each entry pair (p,q) the returned vector contains p:q.
Therefore c(a, b, c) and c(x, y, z) produces c(a:x, b:y, c:z).
}
\usage{multSeq(beg, end)}
\arguments{
    \item{beg}{\code{numeric}. Vector of start positions.}
    %
    \item{end}{\code{numeric}. Vector of end positions.}
}
\value{Integer vector}
\author{Wolfgang Kaisers}
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
% Examples
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
\examples{
# MultSeq
multSeq(c(1, 4, 7), c(2, 5, 8))
}
\keyword{multSeq}
\keyword{segmentize}
