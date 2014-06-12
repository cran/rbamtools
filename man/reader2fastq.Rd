\name{reader2fastq}
\alias{reader2fastq}
\alias{reader2fastq-methods}
\alias{reader2fastq,bamReader-method}
\title{reader2fastq}
\description{The function \code{reader2fastq} writes (remaining) aligns or some part of the stored aligns in a BAM file as reads into a compressed fastq file.}
\usage{reader2fastq(object,filename,which,append=FALSE)}
\arguments{
  \item{object}{bamReader. Must be opened.}
  \item{filename}{Name of output 'fastq' file.}
  \item{which}{(Optional) Logical vector. When given, the routine checks for each of the given values the value of which. When TRUE the next retrieved align will be written to the output file, otherwise will be skipped. }
  \item{append}{(Optional) When TRUE, the routine will append to an existing file. Otherwise existing files will be overwritten.}
}
\value{numeric}
\author{Wolfgang Kaisers}
\examples{
bam<-system.file("extdata","accepted_hits.bam",package="rbamtools")
reader<-bamReader(bam)
reader2fastq(reader,"out.fq.gz")
bamClose(reader)
# Reopen in order to point to first align
reader<-bamReader(bam)
index<-sample(1:100,20)
reader2fastq(reader,"out_subset.fq.gz",which=index)
}
\keyword{reader2fastq}
\keyword{reader2fastq}