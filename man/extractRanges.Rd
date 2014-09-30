\name{extractRanges}
\alias{extractRanges}
\alias{extractRanges-methods}
\alias{extractRanges,bamReader-method}
\title{extractRanges: Extraction of alignments}
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
% Description
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %}
\description{The function \code{extractRanges} takes an opened and indexed
    bamReader object, a \code{data.frame} containg range data, and an ouptut
    filename. 
    From the opened BAM file, aligns in the depicted ranges are transfered
    into a temporary BAM-file.
    The new file will be sorted. The sorted file will be indexed.}
\usage{extractRanges(object,ranges,filename,complex=FALSE,header,idxname)}
\arguments{
    \item{object}{bamReader. Must be opened. Index must be initialized.}
    %
    \item{ranges}{data.frame. Must contain columns 'seqid','start','end'. 
        Each line represents a range on the reference genome for which 
        aligns are extracted.}
    %
    \item{filename}{Character. Name of output BAM-file.
        The filename will be modified.
        The modified filename will have suffix '.bam'.}
    %
    \item{complex}{(default=FALSE) Logical. 
        When TRUE, only aligns with nCigar > 1 are copied.
        When FALSE, all aligns in the depicted ranges are copied.}
    %
    \item{header}{(Optional) bamHeader. It's possible to provide a 
        \code{bamHeader} for the new BAM file (in order to control the 
        header entries in the new file). When no header is given, the header
        from the \code{bamReader} will be used.  
        The 'SO' entry in the 'headerLine' (which gives information
        about the sorting status) segment will be altered.
        The reference sequence dictionary (RSD) of the output file will be
        created new inside the function. It is necessary to do that because
        when there are mismatches between the RSD and the contained aligns,
        it will be impossible to create an index for the new file
        (samtools will crash).}
    %
    \item{idxname}{(Default='filename'.bai) character. 
        The name for the index file.}
}
\author{Wolfgang Kaisers}
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
% Examples
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
\examples{
bam<-system.file("extdata","accepted_hits.bam",package="rbamtools")
reader<-bamReader(bam,idx=TRUE)
# Extract data for HLHL17 gene:
seqid<-"chr1"
start<-895967
end<-901099
gene_name<-"HLHL17" # (optional)
ranges<-data.frame(seqid=seqid,start=start,end=end,gene_name=gene_name)
\dontrun{extractRanges(reader,ranges=ranges,filename="new_file.bam")}
}
\keyword{extractRanges}
\seealso{bamReader}
