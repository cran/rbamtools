\name{bamWriter-class}
\Rdversion{1.1}
\docType{class}
\alias{bamWriter-class}
\alias{bamClose,bamWriter-method}
\alias{bamSave}
\alias{bamSave-methods}
\alias{bamSave,bamWriter-method}
\alias{filename,bamWriter-method}
\alias{initialize,bamWriter-method}
\alias{isOpen,bamWriter-method}
\title{Class \code{"bamWriter"}: Representation of a file connection to
    a BAM file for writing access.}
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
% Description
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
\description{
The bamWriter class represents the writing connection to a BAM-file.
Usually, this class is used to create an empty BAM-file and to copy aligns
from another BAM-file into the new file.
}
\section{Objects from the Class}{
    Objects can be created by calls of the form
    \code{writer<-bamWriter(reader,filename)}.
}
%
\section{Slots}{
    \describe{
        \item{\code{filename}:}{Object of class \code{"character"}.}
        %
        \item{\code{writer}:}{Object of class \code{"externalptr"}.
            Points to samtools struct samtile_t. }
    }
}
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
% Methods
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
\section{Methods}{
    \describe{
        \item{bamClose}{\code{signature(object = "bamWriter")}:
            Closes BAM file. }
        %
        \item{bamSave}{\code{signature(object = "bamWriter",
            value = "bamAlign" or "bamRange",refid = "numeric")}:
            Saves bamAlign or bamRange object to BAM-file.
            refid will be overwritten in all written aligns. }
        %
        \item{filename}{\code{signature(object = "bamWriter")}: 
            Returns filename of opened BAM-file. }
        %
        \item{initialize}{\code{signature(.Object = "bamWriter")}:
            Opens BAM file for writing.}
        %
        \item{isOpen}{\code{signature(object = "bamWriter")}: 
            Returns TRUE when file connection is open. }
    }
}
\references{
\verb{The SAM Format Specification (v1.4-r985). 
The SAM Format Specification Working Group. September 7, 2011.}
\url{http://www.samtools.sourceforge.net/SAM1.pdf}
}
\author{Wolfgang Kaisers}
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
% Examples
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
\examples{
# +++++++++++++++++++++++++++++++++++++++++++++++
# In this example, we copy some complex (i.e. interesting) aligns
# into a new BAM file
bam<-system.file("extdata", "accepted_hits.bam", package="rbamtools")
idx<-paste(bam,"bai",sep=".")

# +++++++++++++++++++++++++++++++++++++++++++++++
# Open BAM file and read complex aligns from chr1
reader<-bamReader(bam)
loadIndex(reader,idx)
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
addSeq(dict,SN="chr10",LN=135534747)
dict
prog<-new("headerProgram")
setVal(prog,"ID","1")
setVal(prog,"PN","tophat")
setVal(prog,"CL","tophat -p8 --library-type fr-unstranded hs_ucsc test.fq")
setVal(prog,"VN","2.0.0")
bh<-bamHeaderText(head=headl,dict=dict,prog=prog)
#getHeaderText(bh)
header<-bamHeader(bh)
# +++++++++++++++++++++++++++++++++++++++++++++++

# +++++++++++++++++++++++++++++++++++++++++++++++
#  Copy aligns in range into new BAM file
\dontrun{
writer<-bamWriter(header, "chr1_complex.bam")
bamSave(writer, range, refid=0)
bamClose(writer)
}
# +++++++++++++++++++++++++++++++++++++++++++++++
}
\keyword{classes}
