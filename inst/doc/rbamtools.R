### R code from vignette source 'rbamtools.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: rbamtools.Rnw:339-343
###################################################
library(rbamtools)
bam<-system.file("extdata","accepted_hits.bam",package="rbamtools")
# Open bam file
reader<-bamReader(bam)


###################################################
### code chunk number 2: rbamtools.Rnw:349-350 (eval = FALSE)
###################################################
## bamSort(reader,prefix="my_sorted",byName=FALSE,maxmem=1e+9)


###################################################
### code chunk number 3: rbamtools.Rnw:354-355 (eval = FALSE)
###################################################
## create.index(reader,idx_filename="index_file_name.bai")


###################################################
### code chunk number 4: rbamtools.Rnw:358-359 (eval = FALSE)
###################################################
## create.index(reader)


###################################################
### code chunk number 5: rbamtools.Rnw:363-365
###################################################
idx<- system.file("extdata", "accepted_hits.bam.bai", package="rbamtools")
load.index(reader,idx)


###################################################
### code chunk number 6: rbamtools.Rnw:368-369
###################################################
index.initialized(reader)


###################################################
### code chunk number 7: rbamtools.Rnw:372-373
###################################################
reader<-bamReader(bam,idx=TRUE)


###################################################
### code chunk number 8: rbamtools.Rnw:379-380
###################################################
getRefData(reader)


###################################################
### code chunk number 9: rbamtools.Rnw:388-392 (eval = FALSE)
###################################################
## header<-getHeader(reader)
## writer<-bamWriter(header,"test.bam")
## # Write aligns using bamSave
## bamClose(writer)


###################################################
### code chunk number 10: rbamtools.Rnw:407-409
###################################################
header<-getHeader(reader)
htxt<-getHeaderText(header)


###################################################
### code chunk number 11: rbamtools.Rnw:428-446
###################################################
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
setVal(prog,"CL","tophat -p8 --library-type fr-unstranded hs_ucsc_index reads.fq")
setVal(prog,"VN","2.0.0")
bh<-bamHeaderText(head=headl,dict=dict,prog=prog)
#getHeaderText(bh)
header<-bamHeader(bh)


###################################################
### code chunk number 12: rbamtools.Rnw:453-454
###################################################
align<-getNextAlign(reader)


###################################################
### code chunk number 13: rbamtools.Rnw:480-491 (eval = FALSE)
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
### code chunk number 14: rbamtools.Rnw:518-529 (eval = FALSE)
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
### code chunk number 15: rbamtools.Rnw:532-533
###################################################
unmapped(align)<-TRUE


###################################################
### code chunk number 16: rbamtools.Rnw:550-554
###################################################
coords<-c(0,899000,900000)
names(coords)<-c("refid","start","stop")
range<-bamRange(reader,coords)
size(range)


###################################################
### code chunk number 17: rbamtools.Rnw:557-562
###################################################
getRefData(reader)
coords<-c(0,0,249250621)
names(coords)<-c("refid","start","stop")
range<-bamRange(reader,coords)
size(range)


###################################################
### code chunk number 18: rbamtools.Rnw:565-569
###################################################
coords<-getRefCoords(reader,"chr1")
coords
range<-bamRange(reader,coords)
size(range)


###################################################
### code chunk number 19: rbamtools.Rnw:578-579
###################################################
align<-getNextAlign(range)


###################################################
### code chunk number 20: rbamtools.Rnw:584-590 (eval = FALSE)
###################################################
## rewind(range)
## while(!is.null(align))
## {
##   # Do whatever is 
##   align<-getNextAlign(range)
## }


###################################################
### code chunk number 21: rbamtools.Rnw:594-595
###################################################
rdf<-as.data.frame(range)


###################################################
### code chunk number 22: rbamtools.Rnw:603-608
###################################################
coords<-getRefCoords(reader,"chr1")
gl<-gapList(reader,coords)
gl
dfr<-as.data.frame(gl)
dfr[1:6,c(1:3,5:8)]


###################################################
### code chunk number 23: rbamtools.Rnw:613-616 (eval = FALSE)
###################################################
## size(gl)
## nAligns(gl)
## nGapAligns(gl)


###################################################
### code chunk number 24: rbamtools.Rnw:649-657
###################################################
coords<-getRefCoords(reader,"chr1")
sl<-siteList(reader,coords)
size(sl)
nAligns(sl)
nGapAligns(sl)
sl
df<-as.data.frame(sl)
head(df)


###################################################
### code chunk number 25: rbamtools.Rnw:667-675
###################################################
bsl<-bamGapList(reader)
bsl
size(bsl)
nAligns(bsl)
nGapAligns(bsl)
summary(bsl)
dfr<-as.data.frame(bsl)
head(dfr)


