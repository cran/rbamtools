### R code from vignette source 'rbamtools.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: rbamtools.Rnw:61-64
###################################################
library(rbamtools)
library(xtable)
options(width=60)


###################################################
### code chunk number 2: rbamtools.Rnw:573-577
###################################################
bam <- system.file("extdata", 
                "accepted_hits.bam", package="rbamtools")
# Open bam file
reader <- bamReader(bam)


###################################################
### code chunk number 3: rbamtools.Rnw:588-590 (eval = FALSE)
###################################################
## bamSort(reader, prefix="my_sorted", 
##             byName=FALSE, maxmem=1e+9)


###################################################
### code chunk number 4: rbamtools.Rnw:597-598 (eval = FALSE)
###################################################
## createIndex(reader, idx_filename="index_file_name.bai")


###################################################
### code chunk number 5: rbamtools.Rnw:605-606 (eval = FALSE)
###################################################
## createIndex(reader)


###################################################
### code chunk number 6: rbamtools.Rnw:613-615
###################################################
idx <- system.file("extdata", "accepted_hits.bam.bai", package="rbamtools")
loadIndex(reader, idx)


###################################################
### code chunk number 7: rbamtools.Rnw:618-619
###################################################
indexInitialized(reader)


###################################################
### code chunk number 8: rbamtools.Rnw:623-624
###################################################
reader <- bamReader(bam, idx=TRUE)


###################################################
### code chunk number 9: rbamtools.Rnw:632-633
###################################################
getRefData(reader)


###################################################
### code chunk number 10: rbamtools.Rnw:646-650 (eval = FALSE)
###################################################
## header <- getHeader(reader)
## writer <- bamWriter(header,"test.bam")
## # Write aligns using bamSave
## bamClose(writer)


###################################################
### code chunk number 11: rbamtools.Rnw:677-679
###################################################
header <- getHeader(reader)
htxt <- getHeaderText(header)


###################################################
### code chunk number 12: rbamtools.Rnw:702-721
###################################################
bh <- new("bamHeaderText")

headl <- new("headerLine")
setVal(headl, "SO", "coordinate")

dict <- new("refSeqDict")
addSeq(dict, SN="chr1",  LN=249250621)
addSeq(dict, SN="chr16", LN=90354753)
dict

prog <- new("headerProgram")
setVal(prog, "ID", "TopHat")
setVal(prog, "PN", "tophat")
setVal(prog, "CL",
    "tophat --library-type fr-unstranded hs_ucsc_index reads.fastq")
setVal(prog, "DS", "Description")
setVal(prog, "VN", "2.0.0")
bh <- bamHeaderText(head=headl, dict=dict, prog=prog)
header <- bamHeader(bh)


###################################################
### code chunk number 13: rbamtools.Rnw:729-730
###################################################
align <- getNextAlign(reader)


###################################################
### code chunk number 14: rbamtools.Rnw:759-770 (eval = FALSE)
###################################################
## name(align)
## flag(align)
## refID(align)
## position(align)
## mapQuality(align)
## cigarData(align)
## nCigar(align)
## mateRefID(align)
## matePosition(align)
## alignSeq(align)
## alignQual(align)


###################################################
### code chunk number 15: rbamtools.Rnw:797-808 (eval = FALSE)
###################################################
## paired(align)
## properPair(align)
## unmapped(align)
## mateUnmapped(align)
## reverseStrand(align)
## mateReverseStrand(align)
## firstInPair(align)
## secondInPair(align)
## secondaryAlign(align)
## failedQC(align)
## pcrORopt_duplicate(align)


###################################################
### code chunk number 16: rbamtools.Rnw:813-814
###################################################
unmapped(align) <- TRUE


###################################################
### code chunk number 17: rbamtools.Rnw:823-832
###################################################
align <- bamAlign("HWUSI-0001", "ATGTACGTCG", "Qual/Strng",
                "4M10N6M", refid=0, position=100)
align
name(align)
alignSeq(align)
alignQual(align)
cigarData(align)
refID(align)
position(align)


###################################################
### code chunk number 18: rbamtools.Rnw:867-871
###################################################
coords <- c(0,899000,900000)
names(coords) <- c("refid","start","stop")
range <- bamRange(reader,coords)
size(range)


###################################################
### code chunk number 19: rbamtools.Rnw:878-883
###################################################
getRefData(reader)
coords <- c(0,0,249250621)
names(coords) <- c("refid","start","stop")
range <- bamRange(reader,coords)
size(range)


###################################################
### code chunk number 20: rbamtools.Rnw:888-892
###################################################
coords <- getRefCoords(reader,"chr1")
coords
range <- bamRange(reader,coords)
size(range)


###################################################
### code chunk number 21: rbamtools.Rnw:902-907
###################################################
range
getCoords(range)
getSeqLen(range)
getParams(range)
getRefName(range)


###################################################
### code chunk number 22: rbamtools.Rnw:913-914
###################################################
getAlignRange(range)


###################################################
### code chunk number 23: rbamtools.Rnw:932-933
###################################################
align <- getNextAlign(range)


###################################################
### code chunk number 24: rbamtools.Rnw:940-946 (eval = FALSE)
###################################################
## rewind(range)
## while(!is.null(align))
## {
##   # Process align data here
##   align <- getNextAlign(range)
## }


###################################################
### code chunk number 25: rbamtools.Rnw:952-953
###################################################
rdf <- as.data.frame(range)


###################################################
### code chunk number 26: rbamtools.Rnw:974-979
###################################################
coords <- getRefCoords(reader, "chr1")
gl <- gapList(reader, coords)
gl
dfr <- as.data.frame(gl)
dfr[1:6, c(1:3, 5:8)]


###################################################
### code chunk number 27: rbamtools.Rnw:989-992 (eval = FALSE)
###################################################
## size(gl)
## nAligns(gl)
## nAlignGaps(gl)


###################################################
### code chunk number 28: rbamtools.Rnw:1036-1044
###################################################
coords <- getRefCoords(reader, "chr1")
sl <- siteList(reader, coords)
size(sl)
nAligns(sl)
nAlignGaps(sl)
sl
df <- as.data.frame(sl)
head(df)


###################################################
### code chunk number 29: rbamtools.Rnw:1059-1067
###################################################
bsl <- bamGapList(reader)
bsl
size(bsl)
nAligns(bsl)
nAlignGaps(bsl)
summary(bsl)
dfr <- as.data.frame(bsl)
head(dfr)


###################################################
### code chunk number 30: rbamtools.Rnw:1093-1096
###################################################
coords <- c(0, 0, 14730)
count <- bamCount(reader, coords)
xtable(matrix(count, nrow=1))


###################################################
### code chunk number 31: rbamtools.Rnw:1100-1101
###################################################
count <- bamCountAll(reader, verbose=TRUE)


###################################################
### code chunk number 32: rbamtools.Rnw:1104-1105
###################################################
xtable(count, digits=0)


###################################################
### code chunk number 33: rbamtools.Rnw:1117-1120
###################################################
align <- bamAlign("HWUSI-0001", "ACCGGGTTTT","Qual/Strng",
                            "4M10N6M", refid=0, position=100)
countNucs(align)


###################################################
### code chunk number 34: rbamtools.Rnw:1123-1127
###################################################
reader <- bamReader(bam, idx=TRUE)
coords <- c(0, 0, 14730)
range <- bamRange(reader, coords)
countNucs(range)


###################################################
### code chunk number 35: rbamtools.Rnw:1143-1144
###################################################
ncs <- nucStats(reader)


###################################################
### code chunk number 36: rbamtools.Rnw:1147-1148
###################################################
xtable(ncs, digits=c(0, 0, 0, 0, 0, 0, 0, 2, 2))


###################################################
### code chunk number 37: rbamtools.Rnw:1156-1157
###################################################
ncs <- nucStats(bam)


###################################################
### code chunk number 38: rbamtools.Rnw:1160-1161
###################################################
xtable(ncs, digits=c(0, 0, 0, 0, 0, 0, 0, 2, 2))


###################################################
### code chunk number 39: rbamtools.Rnw:1182-1183 (eval = FALSE)
###################################################
## createIdxBatch(bam)


###################################################
### code chunk number 40: rbamtools.Rnw:1209-1216 (eval = FALSE)
###################################################
## reader <- bamReader(bam)
## reader2fastq(reader, "out.fastq")
## bamClose(reader)
## # Reopen in order to point to first align
## reader <- bamReader(bam)
## index <- sample(1:100, 20)
## reader2fastq(reader, "out_subset.fastq", which=index)


###################################################
### code chunk number 41: rbamtools.Rnw:1226-1232 (eval = FALSE)
###################################################
## reader <- bamReader(bam,idx=TRUE)
## coords <- as.integer(c(0,0,249250621))
## range <- bamRange(reader,coords)
## range2fastq(range,"rg.fq.gz")
## index <- sample(1:size(range),100)
## range2fastq(range,"rg_subset.fq.gz",which=index)


###################################################
### code chunk number 42: rbamtools.Rnw:1248-1253
###################################################
qdf <- getQualDf(range)
qdf[32:38,1:10]
qdr <- getQualDf(range,prob=TRUE)
qrr <- round(qdr,2)
qrr[32:38,1:10]


###################################################
### code chunk number 43: rbamtools.Rnw:1260-1262
###################################################
qt <- getQualQuantiles(range,c(0.25,0.5,0.75))
qt[,1:10]


###################################################
### code chunk number 44: rbamtools.Rnw:1268-1269
###################################################
plotQualQuant(range)


###################################################
### code chunk number 45: rbamtools.Rnw:1288-1305
###################################################
# WASH7P coordinates
xlim <- c(10000, 30000)
coords <- c(0,xlim[1], xlim[2])
range <- bamRange(reader, coords)
bamClose(reader)
ad <- alignDepth(range)
ad
getParams(ad)
# Identifier
gene <- "WASH7P"
ensg_id <- "ENSG00000227232"
enst_id <- "ENST00000538476"
# Get exon positions
start <- c(14411, 15000, 15796, 15904, 16607, 16748, 16858, 17233,
                    17602, 17915, 18268, 24737, 29534)
end <-   c(14502, 15038, 15901, 15947, 16745, 16765, 17055, 17364,
                    17742, 18061, 18366, 24891, 29806)


###################################################
### code chunk number 46: rbamtools.Rnw:1308-1316
###################################################
plotAlignDepth(ad, lwd = 2, xlim = xlim,
            main = paste("Align depth for gene",gene),
            ylab = "Align depth", start = start,
            end = end, strand = "-",
            transcript = paste("Chromosome 1",
                "\tGene ENSG00000227232", ensg_id, 
                "\tTranscript ",enst_id
))


