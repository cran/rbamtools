\name{getVal}
\alias{getVal}
\title{getVal(object,members)}
\description{Retrieving value of data-member from object}
\usage{getVal(object,member)}
\arguments{
  \item{object}{An instance of class headerLine, headerReadGroup or headerProgram }
  \item{member}{A (character) member identifier corresponding to the data-members of the given object (e.g. PN for headerProgram)}
}
\details{
Data members for class headerLine:
\tabular{ll}{
VN \tab Format version \cr
SO \tab Sorting order  \cr
}
Valid values for sorting order (SO) are: unknown (default), unsorted, queryname or coordinate.

Data members for class readGroup:
\tabular{ll}{
ID \tab Read Group identifier \cr
CN \tab Name of sequencing center \cr
DS \tab Description \cr
FO \tab Flow order \cr
KS \tab Nucleotides corresponding to key sequence of each read \cr
LB \tab Library \cr
PG \tab Programs used for processing the Read Group \cr
PI \tab Predicted median insert size \cr
PL \tab Sequencing Platform \cr
SM \tab Sample name \cr
}
Valid values for Sequencing Platform (PL) are: CAPILLARY, LS454, ILLUMINA, SOLID, HELICOS, IONTORRENT or PACBIO.

Data members for class headerProgram
\tabular{ll}{
ID \tab Program record identifier \cr
PN \tab Program name \cr
CL \tab Command line \cr
PP \tab Previous @PG-ID (Must match another ID in list) \cr
VN \tab Program Version \cr
}
Data members for refSeqDict can be accessed via data.frame generic functions (head,tail,[,[<-) or by conversion (as.data.frame).
}
\author{Wolfgang Kaisers}
\examples{
library(rbamtools)
bam<-system.file("extdata","accepted_hits.bam",package="rbamtools")
reader<-bamReader(bam)
header<-getHeader(reader)
htxt<-getHeaderText(header)
hl<-headerLine(htxt)
getVal(hl,"SO")
bamClose(reader)
}
\keyword{getVal}
