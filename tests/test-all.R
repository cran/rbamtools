

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# Load prerequisites
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #

require(rbamtools)

##============================================================================##
## A    Initialize example data
##============================================================================##

bam<- system.file("extdata", "accepted_hits.bam", package="rbamtools")
idx<- system.file("extdata", "accepted_hits.bam.bai", package="rbamtools")
# Open BAM-file for reading
reader<-bamReader(bam,idx=TRUE)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# Run tests
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #

##============================================================================##
## B    Bam align
##============================================================================##
align <- getNextAlign(reader)


##----------------------------------------------------------------------------##
## B.1 Test for accessing new align FLAG (0x800)
##----------------------------------------------------------------------------##

f <- flag(align)
suppAlign(align)
suppAlign(align) <- TRUE
if(!suppAlign(align))
    stop("[test_bam_align.r] Setting suppAlign failed")

suppAlign(align) <- FALSE
if(flag(align)!=f)
    stop("[test_bam_align.r]
            Setting and re-setting suppAlign changed other flags")


##============================================================================##
## B    Bam header
##============================================================================##

header<-getHeader(reader)


##----------------------------------------------------------------------------##
## B.1 headerReadGroup
##----------------------------------------------------------------------------##

##                                                                            ##
## B.1.1    Create headerReadGroup by simply parsing text
##                                                                            ##

hrg <- c("@RG\tID:rg1\tCN:seqCenter1\tDT:01.01.2011\tSM:sm1",
         "@RG\tID:rg2\tCN:seqCenter2\tDT:01.01.2012\tSM:sm2")
object <- new("headerReadGroup", hrg)

if(!identical(object@ID,paste("rg", 1:2, sep="")))
    stop("[test_bam_header.r] Setting ID slot in headerReadGroup failed")


##                                                                            ##
## B.1.2 Create headerReadGroup by user interface  and check for consistency
##          of returned header text
##                                                                            ##

hrg <- new("headerReadGroup")
hrg <- addReadGroup(hrg, list(ID="rg1", CN="sct1", FO="fo1"))
hrg <- addReadGroup(hrg, list(ID="rg2", CN="sct2", FO="fo2", LB="lb2"))
hrg <- addReadGroup(hrg, list(ID="rg3", CN="sct3", LB="lb3"))
hrg2 <- new("headerReadGroup", getHeaderText(hrg))

if(!identical(hrg, hrg2))
    stop("[test_bam_header.r] Adding headerReadGroup's failed")

##----------------------------------------------------------------------------##
## B.2 headerProgram
##----------------------------------------------------------------------------##

prog <- new("headerProgram")
setVal(prog, "ID", "pg.001")
setVal(prog, "PN", "Program name")
setVal(prog, "CL", "Command line")
setVal(prog, "DS", "Description")

if(!identical(getVal(prog,"ID"), "pg.001"))
    stop("[test_bam_header.r] Setting ID in headerProgram failed")

if(!identical(getVal(prog,"DS"), "Description"))
    stop("[test_bam_header.r] Setting Description in header Program failed")

##                                                                            ##
## B.2.1 Set new Description segment
##                                                                            ##
header<-getHeader(reader)
htxt<-getHeaderText(header)
prog <-headerProgram(htxt)

setVal(prog, "DS", "Description")
headerProgram(htxt) <- prog

##                                                                            ##
## B.2.2    Convert to binary header representation and retrieve Description
##          segment
##                                                                            ##
bh <- bamHeader(htxt)
if(!identical(getVal(headerProgram(getHeaderText(bh)),"DS"),"Description"))
    stop("[test_bam_header.r] Writing Program description to bamHeader failed")



##============================================================================##
## C    BamRange
##============================================================================##
coords<-getRefCoords(reader,"chr1")
rg<-bamRange(reader,coords)


if(size(rg)!=2216)
    stop("[test_bam_range.r] size of range must be 2216!")



##============================================================================##
## D    Genome-Partition
##============================================================================##

##----------------------------------------------------------------------------##
## D.1  Create genomePartition object
##----------------------------------------------------------------------------##
# Open (indexed) BAM file (done by calling test-all.R)
# bam<-system.file("extdata", "accepted_hits.bam", package="rbamtools")
# reader<-bamReader(bam,idx=TRUE)

# Provide exon positions
id <- 1:13
seqid <- "chr1"
gene <- "WASH7P"
ensg_id <- "ENSG00000227232"
start <- c(14411, 15000, 15796, 15904, 16607, 16748, 16858, 17233,
           17602, 17915, 18268, 24737, 29534)
end <-   c(14502, 15038, 15901, 15947, 16745, 16765, 17055, 17364,
           17742, 18061, 18366, 24891, 29806)

ref <- data.frame(id=id, seqid=seqid, begin=start, end=end, gene=gene, ensg=ensg_id)

##----------------------------------------------------------------------------##
## D.2  Create partition (adds equidistant grid)
##----------------------------------------------------------------------------##
partition <- genomePartition(reader, ref)

# data_frame test:
rn <- as.numeric(rownames(partition@ev$reflist[[1]]))
idx <- 1:length(rn)
if(any(rn!=idx))
    stop("[test_genome_partition.r] Error in rownames produced by data_frame!")



##============================================================================##
## E    Count alignments in specified genomic segments
##============================================================================##
coords <- c(0, 0, 2e4)
segments <- seq(14000, 20000, 20)
segcount<-rangeSegCount(reader, coords, segments)
segcount
dfr<-as.data.frame(segcount)

if(nrow(dfr)!=301)
    stop("[test_range_seg_count.r] data.frame has wrong number of lines!")

if(sum(dfr$count)!=2112)
    stop("[test_range_seg_count.r] Wrong total number of aligns!")

##============================================================================##
## F    Test vector segmentation
##============================================================================##

if(any(multSeq(c(1, 4, 7), c(2, 5, 8)) != c(1,2, 4,5, 7,8)))
    stop("[test_vector_segmentation] Wrong multSeq value!")

if(any(segmentize(1:10, c(1, 4, 7), c(2, 5, 8)) != c(1,2, 4,5, 7,8)))
    stop("[test_vector_segmentation] Wrong segmentize.ANY value!")

if(any(segmentize(letters, c(1,4,7), c(2,5,8)) != c("a","b","d","e","g","h")))
    stop("[test_vector_segmentation] Wrong segmentize.letters value!")

# Create named input vector ()
x <- rep(0, 11)
names(x) <- 10:20
x[3:5] <- 1:3
x[7:9] <- 4:6

if(any(segmentize(x, c(12, 16), c(14, 18), offset=10)!=1:6))
    stop("[test_vector_segmentation] Wrong output value!")


##============================================================================##
## G    Cleanup
##============================================================================##

bamClose(reader)
rm(reader)
gc()

cat("[rbamtools] rest-all.R tests finished.\n")


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# END OF FILE
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
