\name{rbamtools-package}
\alias{rbamtools-package}
\alias{rbamtools}
\docType{package}
\title{Reading, writing and manipulating BAM-file format.}
\description{
BAM is Binary (Sequence) Alignment/Map format. Many sequence alignment products which align Second Generation Sequence reads to a Sequence reference (such as the human genome) use BAM-file format as output.
Analysis of results of a sequence alignment requires reading and interpreting BAM-files and sometimes manipulating BAM-files. The rbamtools package provides a R-interface to the samtools C-library by Heng Li.
}
\details{
\tabular{ll}{
Package: \tab rbamtools\cr
Type: \tab Package\cr
Version: \tab 2.0\cr
Date: \tab 2012-02-03\cr
License: \tab Artistic 2.0 + MIT License (see LICENSE in src/samtools subdirectory)\cr
LazyLoad: \tab yes\cr
Depends: \tab methods\cr
}

The package is organized in S4 classes. Four classes represent the data and organize the defined functionality:
\tabular{ll}{
bamReader \tab Reading aligns from a BAM-file\cr
bamAlign  \tab Single BAM alignment (including Data accessor functions)\cr
bamRange  \tab (Linked-) List container for BAM alignments\cr
bamWriter \tab Writing alings to a BAM-file\cr
gapList   \tab (Linked-) List container for alignment gaps.\cr
}
}
\author{
Wolfgang Kaisers
Maintainer: Wolfgang Kaisers <kaisers@med.uni-duesseldorf.de>
}
\references{
\verb{The SAM Format Specification (v1.4-r985). The SAM Format Specification Working Group. September 7, 2011.}
\url{http://www.samtools.sourceforge.net/SAM1.pdf}

\verb{
Li H.*, Handsaker B.*, Wysoker A., Fennell T., Ruan J., Homer N., Marth G., Abecasis G., Durbin R. and 1000 Genome Project Data Processing Subgroup (2009)
The Sequence alignment/map (SAM) format and SAMtools.
Bioinformatics, 25, 2078-9.
}
}

\keyword{package}
\examples{
library(rbamtools)
bam <- system.file("extdata", "accepted_hits.bam", package="rbamtools")
idx <- system.file("extdata", "accepted_hits.bam.bai", package="rbamtools")
reader<-bamReader(bam)
isOpen(reader)
align<-getNextAlign(reader)
load.index(reader,idx)
index.initialized(reader)
coords<-as.integer(c(0,0,249250621))
range<-bamRange(reader,coords)
align<-getNextAlign(range)
bamClose(reader)
}
