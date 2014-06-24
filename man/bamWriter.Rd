\name{bamWriter}
\alias{bamWriter}
\title{bamWriter Object}
\description{The bamWriter function takes a \code{bamHeader} object and a filename and returns an object of class 'bamWriter' which represents a writing connection to a BAM-file}
\usage{
bamWriter(x,filename)
}
\arguments{
  \item{x}{An instance of class bamHeader.}
  \item{filename}{Filename of BAM-file to be opened for writing.}
}
\author{Wolfgang Kaisers}
\examples{
# +++++++++++++++++++++++++++++++++++++++++++++++
# In this example, we copy some complex (i.e. interesting) aligns
# into a new BAM file
bam<-system.file("extdata", "accepted_hits.bam", package="rbamtools")
idx<-paste(bam,"bai",sep=".")

# +++++++++++++++++++++++++++++++++++++++++++++++
# Open BAM file and read complex aligns from chr1
reader<-bamReader(bam)
load.index(reader,idx)
coords<-as.integer(c(0,0,249250621))
range<-bamRange(reader,coords,complex=TRUE)
bamClose(reader)
# +++++++++++++++++++++++++++++++++++++++++++++++

# +++++++++++++++++++++++++++++++++++++++++++++++
#  Create bamHeader from scratch
bh<-new("bamHeaderText")
headl<-new("headerLine")
setVal(headl,"SO","coordinate")
dict<-new("refSeqDict")
addSeq(dict,SN="chr1",LN=249250621)
dict
prog<-new("headerProgram")
setVal(prog,"ID","1")
setVal(prog,"PN","tophat")
setVal(prog,"CL","tophat -p8 --library-type fr-unstranded hs_ucsc rna033.fastq")
setVal(prog,"VN","2.0.0")
bh<-bamHeaderText(head=headl,dict=dict,prog=prog)
#getHeaderText(bh)
header<-bamHeader(bh)
# +++++++++++++++++++++++++++++++++++++++++++++++

# +++++++++++++++++++++++++++++++++++++++++++++++
#  Copy aligns in range into new BAM file
writer<-bamWriter(header,"chr1_complex.bam")
bamSave(writer,range,refid=0)
bamClose(writer)
# +++++++++++++++++++++++++++++++++++++++++++++++
}
\keyword{bamWriter}
\keyword{bamAlign}
