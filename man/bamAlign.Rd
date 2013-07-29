\name{bamAlign}
\alias{bamAlign}
\title{bamAlign}
\description{The function takes data stored in aling-fields and creates a \code{bamAlign} object, which can be stored in a BAM file via bamWriter.}
\usage{bamAlign(qname,qseq,qqual,cigar,refid,position,flag=272L,alqual=10L,
                   mrefid=(-1L),mpos=(-1L),insertsize=0L)}
\arguments{
  \item{qname}{Query name, e.g. "HWUSI..." for Illumina sequences.}
  \item{qseq}{Query sequence (DNA-sequence)}
  \item{qqual}{Query quality (ASCII coded quality values). Must contain same number of characters as qseq}
  \item{cigar}{CIGAR string. Must be in valid format, e.g. 45M100N56N. Sequence length must match encoded items in CIGAR string.}
  \item{refid}{Integer. 0-based index which must have a counterpart in Reference Sequence Dictionary (otherwise samtools crashed when creating a BAM-index file).}
  \item{position}{Integer.}
  \item{flag}{Integer. Contains information about binary stored flags in align.}
  \item{alqual}{Integer. Alignment quality.}
  \item{mrefid}{Integer. Mate refid. Used for paired end reads.}
  \item{mpos}{Integer. Mate position. Used for paired end reads.}
  \item{insertsize}{Integer.}
}
\value{bamGapList}
\author{Wolfgang Kaisers}
\examples{
align<-bamAlign("HWUSI-0001","ATGTACGTCG","Qual/Strng","4M10N6M",refid=0,position=100)
align
name(align)
alignSeq(align)
alignQual(align)
cigarData(align)
refID(align)
position(align)
}
\keyword{bamAlign}