### R code from vignette source 'rbamtools.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: rbamtools.Rnw:365-369
###################################################
library(rbamtools)
bam<-system.file("extdata","accepted_hits.bam",package="rbamtools")
# Open bam file
reader<-bamReader(bam)


###################################################
### code chunk number 2: rbamtools.Rnw:375-376 (eval = FALSE)
###################################################
## bamSort(reader,prefix="my_sorted",byName=FALSE,maxmem=1e+9)


###################################################
### code chunk number 3: rbamtools.Rnw:380-381 (eval = FALSE)
###################################################
## create.index(reader,idx_filename="index_file_name.bai")


###################################################
### code chunk number 4: rbamtools.Rnw:384-385 (eval = FALSE)
###################################################
## create.index(reader)


###################################################
### code chunk number 5: rbamtools.Rnw:389-391
###################################################
idx<- system.file("extdata", "accepted_hits.bam.bai", package="rbamtools")
load.index(reader,idx)


###################################################
### code chunk number 6: rbamtools.Rnw:394-395
###################################################
index.initialized(reader)


###################################################
### code chunk number 7: rbamtools.Rnw:398-399
###################################################
reader<-bamReader(bam,idx=TRUE)


###################################################
### code chunk number 8: rbamtools.Rnw:405-406
###################################################
getRefData(reader)


###################################################
### code chunk number 9: rbamtools.Rnw:414-418 (eval = FALSE)
###################################################
## header<-getHeader(reader)
## writer<-bamWriter(header,"test.bam")
## # Write aligns using bamSave
## bamClose(writer)


###################################################
### code chunk number 10: rbamtools.Rnw:433-435
###################################################
header<-getHeader(reader)
htxt<-getHeaderText(header)


###################################################
### code chunk number 11: rbamtools.Rnw:454-472
###################################################
bh<-new("bamHeaderText")

headl<-new("headerLine")
setVal(headl,"SO","coordinate")

dict<-new("refSeqDict")
addSeq(dict,SN="chr1",LN=249250621)
addSeq(dict,SN="chr16",LN=90354753)
dict

prog<-new("headerProgram")
setVal(prog,"ID","1")
setVal(prog,"PN","tophat")
setVal(prog,"CL","tophat --library-type fr-unstranded hs_ucsc_index reads.fastq")
setVal(prog,"VN","2.0.0")
bh<-bamHeaderText(head=headl,dict=dict,prog=prog)
#getHeaderText(bh)
header<-bamHeader(bh)


###################################################
### code chunk number 12: rbamtools.Rnw:479-480
###################################################
align<-getNextAlign(reader)


###################################################
### code chunk number 13: rbamtools.Rnw:506-517 (eval = FALSE)
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
### code chunk number 14: rbamtools.Rnw:544-555 (eval = FALSE)
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
### code chunk number 15: rbamtools.Rnw:558-559
###################################################
unmapped(align)<-TRUE


###################################################
### code chunk number 16: rbamtools.Rnw:564-572
###################################################
align<-bamAlign("HWUSI-0001","ATGTACGTCG","Qual/Strng","4M10N6M",refid=0,position=100)
align
name(align)
alignSeq(align)
alignQual(align)
cigarData(align)
refID(align)
position(align)


###################################################
### code chunk number 17: rbamtools.Rnw:589-593
###################################################
coords<-c(0,899000,900000)
names(coords)<-c("refid","start","stop")
range<-bamRange(reader,coords)
size(range)


###################################################
### code chunk number 18: rbamtools.Rnw:596-601
###################################################
getRefData(reader)
coords<-c(0,0,249250621)
names(coords)<-c("refid","start","stop")
range<-bamRange(reader,coords)
size(range)


###################################################
### code chunk number 19: rbamtools.Rnw:604-608
###################################################
coords<-getRefCoords(reader,"chr1")
coords
range<-bamRange(reader,coords)
size(range)


###################################################
### code chunk number 20: rbamtools.Rnw:612-617
###################################################
range
getCoords(range)
getSeqLen(range)
getParams(range)
getRefName(range)


###################################################
### code chunk number 21: rbamtools.Rnw:620-621
###################################################
getAlignRange(range)


###################################################
### code chunk number 22: rbamtools.Rnw:631-632
###################################################
align<-getNextAlign(range)


###################################################
### code chunk number 23: rbamtools.Rnw:637-643 (eval = FALSE)
###################################################
## rewind(range)
## while(!is.null(align))
## {
##   # Process align data here
##   align<-getNextAlign(range)
## }


###################################################
### code chunk number 24: rbamtools.Rnw:647-648
###################################################
rdf<-as.data.frame(range)


###################################################
### code chunk number 25: rbamtools.Rnw:656-661
###################################################
coords<-getRefCoords(reader,"chr1")
gl<-gapList(reader,coords)
gl
dfr<-as.data.frame(gl)
dfr[1:6,c(1:3,5:8)]


###################################################
### code chunk number 26: rbamtools.Rnw:666-669 (eval = FALSE)
###################################################
## size(gl)
## nAligns(gl)
## nAlignGaps(gl)


###################################################
### code chunk number 27: rbamtools.Rnw:702-710
###################################################
coords<-getRefCoords(reader,"chr1")
sl<-siteList(reader,coords)
size(sl)
nAligns(sl)
nAlignGaps(sl)
sl
df<-as.data.frame(sl)
head(df)


###################################################
### code chunk number 28: rbamtools.Rnw:717-725
###################################################
bsl<-bamGapList(reader)
bsl
size(bsl)
nAligns(bsl)
nAlignGaps(bsl)
summary(bsl)
dfr<-as.data.frame(bsl)
head(dfr)


###################################################
### code chunk number 29: rbamtools.Rnw:736-742
###################################################
bam<-system.file("extdata","accepted_hits.bam",package="rbamtools")
coords<-c(0,0,14730)
count<-bamCount(reader,coords)
count
count<-bamCountAll(reader,verbose=TRUE)
count


###################################################
### code chunk number 30: rbamtools.Rnw:748-750
###################################################
align<-bamAlign("HWUSI-0001","ACCGGGTTTT","Qual/Strng","4M10N6M",refid=0,position=100)
countNucs(align)


###################################################
### code chunk number 31: rbamtools.Rnw:753-758
###################################################
bam<-system.file("extdata","accepted_hits.bam",package="rbamtools")
reader<-bamReader(bam,idx=TRUE)
coords<-c(0,0,14730)
range<-bamRange(reader,coords)
countNucs(range)


###################################################
### code chunk number 32: rbamtools.Rnw:765-766
###################################################
nucStats(reader)


###################################################
### code chunk number 33: rbamtools.Rnw:770-771
###################################################
nucStats(bam)


###################################################
### code chunk number 34: rbamtools.Rnw:780-782 (eval = FALSE)
###################################################
## bam<-system.file("extdata","accepted_hits.bam",package="rbamtools")
## create.idx.batch(bam)


###################################################
### code chunk number 35: rbamtools.Rnw:792-800 (eval = FALSE)
###################################################
## bam<-system.file("extdata","accepted_hits.bam",package="rbamtools")
## reader<-bamReader(bam)
## reader2fastq(reader,"out.fastq")
## bamClose(reader)
## # Reopen in order to point to first align
## reader<-bamReader(bam)
## index<-sample(1:100,20)
## reader2fastq(reader,"out_subset.fastq",which=index)


###################################################
### code chunk number 36: rbamtools.Rnw:805-812 (eval = FALSE)
###################################################
## bam<-system.file("extdata","accepted_hits.bam",package="rbamtools")
## reader<-bamReader(bam,idx=TRUE)
## coords<-as.integer(c(0,0,249250621))
## range<-bamRange(reader,coords)
## range2fastq(range,"rg.fq.gz")
## index<-sample(1:size(range),100)
## range2fastq(range,"rg_subset.fq.gz",which=index)


###################################################
### code chunk number 37: rbamtools.Rnw:821-826
###################################################
qdf<-getQualDf(range)
qdf[32:38,1:10]
qdr<-getQualDf(range,prob=TRUE)
qrr<-round(qdr,2)
qrr[32:38,1:10]


###################################################
### code chunk number 38: rbamtools.Rnw:831-833
###################################################
qt<-getQualQuantiles(range,c(0.25,0.5,0.75))
qt[,1:10]


###################################################
### code chunk number 39: rbamtools.Rnw:838-839
###################################################
plotQualQuant(range)


###################################################
### code chunk number 40: rbamtools.Rnw:851-858
###################################################
# WASH7P coordinates
coords<-as.integer(c(0,16950,17400))
range<-bamRange(reader,coords)
bamClose(reader)
ad<-alignDepth(range)
ad
getParams(ad)


###################################################
### code chunk number 41: rbamtools.Rnw:861-862
###################################################
plotAlignDepth(ad,col="lightblue")


