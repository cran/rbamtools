
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #
#  File   : rbamtools.r                                                         #
#  Date   : 21.Sep.2012                                                         #
#  Content: R-Source for package rbamtools                                      #
#  Version: 2.1.5                                                               #
#  Author : W. Kaisers                                                          #
#  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #
#  CRAN submission:
#  wput rbamtools_2.0.tar.gz ftp://cran.r-project.org/incoming/ 
#  Changelog
#  30.Okt.12  [initialize.gapList] Made printout message optional (verbose)
#  31.Okt.12  [bamRange] Included test for initialized index
#  31.Okt.12  Check for open reader in getHeader, getHeaderText, getRefCount
#  01.Nov.12  [get_const_next_align] added to correct memory leak.
#  08.Nov.12  Reading and writing big bamRanges (pure C, no R) valgrind checked.
#  09.Nov.12  [bamCopy.bamReader] Added which allows refwise copying.
#  31.Dec.12  gapSiteList class added
#  11.Jan.13  bamSiteList class added
#  06.Feb.13  First successful test of bamSiteList on 36 BAM-files (871.926/sec)
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #

#.Last.lib<-function(libpath) { library.dynam.unload("rbamtools",libpath) }
.onUnload<-function(libpath) { library.dynam.unload("rbamtools",libpath) }

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #
#  Declaration of generics
setGeneric("filename", function(object) standardGeneric("filename"))
setGeneric("isOpen",function(con,rw="") standardGeneric("isOpen"))
setGeneric("bamClose", function(object) standardGeneric("bamClose"))
# generic for bamReader and bamRange
setGeneric("getNextAlign",function(object) standardGeneric("getNextAlign")) 
# generic for bamWriter and bamRange
setGeneric("bamSave",function(object,...) standardGeneric("bamSave"))
# Generic for conversion into list
setGeneric("as.list",function(x,...) standardGeneric("as.list"))
# Generic for retrieving RefData string from Objects
setGeneric("getHeaderText",function(object,delim="\n") standardGeneric("getHeaderText"))
# Generic for Reading member from object list
setGeneric("getVal",function(object,member)standardGeneric("getVal"))
# Generic for Writing member to object list
setGeneric("setVal",function(object,members,values)standardGeneric("setVal"))
# Generic for retrieving of list size
setGeneric("size",function(object) standardGeneric("size"))
# Generic for retrieving Nr of aligns in BAM region from gapList
setGeneric("nAligns",   function(object)standardGeneric("nAligns"))
# Generic for retrieving Nr of gapped-aligns in BAM region from gapList
setGeneric("nGapAligns",function(object)standardGeneric("nGapAligns"))
# Generic for reading gapLists (align gaps) from bamReader
setGeneric("gapList",function(object,coords)standardGeneric("gapList"))
# Generic for reading gapSiteList (merged align gap sites) from bamReader
setGeneric("siteList",function(object,coords)standardGeneric("siteList"))
# Generic for reading bamSiteList (merged align gap sites for whole bam-files) from bamReader
setGeneric("bamSiteList",function(object)standardGeneric("bamSiteList"))


###################################################################################################
#                                                                                                 #
# bamReader                                                                                       #
#                                                                                                 #
###################################################################################################

setClass("bamReader",representation(filename="character",reader="externalptr",
				index="externalptr",startpos="numeric"),
         validity=function(object){return(ifelse(is.null(object@reader),FALSE,TRUE))})

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #
#  Opening and closing a BAM-File for reading
#  
setMethod(f="initialize", signature="bamReader",
          definition=function(.Object,filename){
            .Object@filename<-filename
            .Object@reader=.Call("bam_reader_open",path.expand(filename),PACKAGE="rbamtools")
            .Object@startpos=.Call("bam_reader_tell",.Object@reader,PACKAGE="rbamtools")
            return(.Object)
          }
)

bamReader<-function(filename,indexname,idx=FALSE,verbose=0){
  if(!is.logical(idx))
    stop("[bamReader] idx must be logical!")
  if(length(idx)>1)
    stop("[bamReader] length(idx) must be 1!")
  if(!is.numeric(verbose))
    stop("[bamReader] verbose must be numeric!")
    
  reader<-new("bamReader",filename)
  if((!idx) && missing(indexname))
  {
    if(verbose[1]==1)
      cat("[bamReader] Opened file '",basename(filename),"'.\n",sep="")
    else if(verbose[1]==2)
      cat("[bamReader] Opened file '",filename,"'.\n",sep="")
    return(reader)
  }
  
  # use indexname or set default
  if(missing(indexname))
    idxfile<-paste(filename,"bai",sep=".")
  else
    idxfile<-indexname   
  loadIndex(reader,idxfile)
 
  if(verbose[1]==1)
    cat("[bamReader] Opened file '",basename(filename),"' and index '",basename(idxfile),"'.\n",sep="")
  else if(verbose[1]==2)
    cat("[bamReader] Opened file '",filename,"' and index '",idxfile,"'.\n",sep="")
  
  return(reader)
}

setMethod("filename", "bamReader", function(object) return(object@filename))
setMethod("isOpen",signature="bamReader",definition=function(con,rw="")
{ return(!(.Call("is_nil_externalptr",con@reader,PACKAGE="rbamtools"))) })
setMethod(f="bamClose",signature="bamReader",definition=function(object) {
  if(!.Call("is_nil_externalptr",object@index,PACKAGE="rbamtools"))
  {.Call("bam_reader_unload_index",object@index,PACKAGE="rbamtools")}
  invisible(.Call("bam_reader_close",object@reader,PACKAGE="rbamtools"))
})

#  End: Opening and closing a BAM-File for reading
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #
#  Header related functions

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #
#  This is one standard Method for creation of bamHeader                        #
#  and is used as a simple way to pass a header to a new                        #
#  instance of bamWriter                                                        #
setGeneric("getHeader",function(object)standardGeneric("getHeader"))
setMethod(f="getHeader",signature="bamReader",definition=function(object){
  if(!isOpen(object))
    stop("[getHeader.bamReader] reader must be opened! Check with 'isOpen(reader)'!")
  return(new("bamHeader",.Call("bam_reader_get_header",object@reader))) })
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #

setMethod(f="getHeaderText",signature="bamReader",definition=function(object){
  if(!isOpen(object))
    stop("[getHeaderText.bamReader] reader must be opened! Check with 'isOpen(reader)'!")
  return(new("bamHeaderText",
             .Call("bam_reader_get_header_text",object@reader,PACKAGE="rbamtools")))
})

# getRefCount
setGeneric("getRefCount",function(object) standardGeneric("getRefCount"))
setMethod(f="getRefCount",signature="bamReader",definition=function(object) {
  if(!isOpen(object))
    stop("[getRefCount.bamReader] reader must be opened! Check with 'isOpen(reader)'!")
  return(.Call("bam_reader_get_ref_count",object@reader,PACKAGE="rbamtools"))})

# getRefData
setGeneric("getRefData",function(object) standardGeneric("getRefData"))
setMethod(f="getRefData",signature="bamReader",definition=function(object) {
  if(!isOpen(object))
    stop("[getRefData.bamReader] reader must be opened! Check with 'isOpen(reader)'!")
  return(.Call("bam_reader_get_ref_data",object@reader,PACKAGE="rbamtools"))})

# getRefCoords: Helper function that returns coordinates of entire ref
# for usage with bamRange, gapList or siteList function
setGeneric("getRefCoords",function(object,sn)standardGeneric("getRefCoords"))
setMethod(f="getRefCoords",signature="bamReader",definition=function(object,sn){
  if(!is.character(sn))
    stop("[getRefCoords] sn must be character!")
  if(length(sn)>1)
    stop("[getRefCoords] sn must have length 1!")
  ref<-getRefData(object)
  id<-which(sn==ref$SN)
  if(length(id)==0)
    stop("[getRefCoords] No match for sn in ref-data-SN!")
  return(c(ref$ID[id],0,ref$LN[id]))
})


#  End Header related functions
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #
#  Index related functions

# createIndex
setGeneric("createIndex",function(object,idx_filename) standardGeneric("createIndex"))
setMethod(f="createIndex",signature="bamReader",definition=function(object,idx_filename)
{invisible(.Call("bam_reader_create_index",path.expand(object@filename),
                 path.expand(idx_filename),PACKAGE="rbamtools"))})

setGeneric("loadIndex",function(object,filename) standardGeneric("loadIndex"))
setMethod("loadIndex",signature="bamReader",definition=function(object,filename){
  if(!is.character(filename))
    stop("[bamReader.loadIndex] Filename must be character!\n")
  if(!file.exists(filename))
    stop("[bamReader.loadIndex] Index file \"",filename,"\" does not exist!\n")
  
  # Set index Variable in given bamReader object:
  # Read object name, create expression string and evaluate in parent frame
  reader<-deparse(substitute(object))
  extxt<-paste(reader,"@index<-.Call(\"bam_reader_load_index\",\"",
               path.expand(filename),"\",PACKAGE=\"rbamtools\")",sep="")
  eval.parent(parse(text=extxt))
  
  # Return true if bamReader@index!=NULL (parent frame)
  extxt<-paste(".Call(\"is_nil_externalptr\",",reader,"@index,PACKAGE=\"rbamtools\")",sep="")
  invisible(!eval.parent(parse(text=extxt)))
})
setGeneric("index.initialized",function(object) standardGeneric("index.initialized"))
setMethod("index.initialized", signature="bamReader",definition=function(object)
{return(!(.Call("is_nil_externalptr",object@index,PACKAGE="rbamtools")))})

setGeneric("bamSort",function(object,prefix,byName=FALSE,maxmem=1e+9) standardGeneric("bamSort"))
setMethod(f="bamSort",signature="bamReader",
          definition=function(object,prefix="sorted",byName=FALSE,maxmem=1e+9)
          {
            maxmem<-floor(maxmem)
            cat("[bamSort] Filename: ",object@filename,"\n")
            cat("[bamSort] Prefix  : ",prefix,"\n")
            cat("[bamSort] Maxmem  : ",maxmem,"\n")
            cat("[bamSort] By Name : ",byName,"\n")
            .Call("bam_reader_sort_file",object@filename,prefix,maxmem,byName,PACKAGE="rbamtools")
            cat("[bamSort] Sorting finished.\n")
            return(invisible(paste(prefix,".bam",sep="")))
          })

#  End Index related functions
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #


# getNextAlign
setMethod(f="getNextAlign",signature="bamReader",definition=function(object)
{
  ans<-.Call("bam_reader_get_next_align",object@reader,PACKAGE="rbamtools")
  if(is.null(ans))
    return(invisible(NULL))
  else
    return(new("bamAlign",ans))
})


# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #
# Reading gap-lists
setMethod("gapList","bamReader",function(object,coords)
{
  if(!index.initialized(object))
    stop("[gapList.bamReader] Reader must have initialized index!")
  return(new("gapList",object,coords))
})
setMethod("siteList","bamReader",function(object,coords)
{
  if(!index.initialized(object))
    stop("[siteList.bamReader] Reader must have initialized index!")
  return(new("gapSiteList",object,coords))
})
setMethod("bamSiteList","bamReader",function(object)
{
  if(!index.initialized(object))
    stop("[bamSiteList.bamReader] Reader must have initialized index!")
  return(new("bamSiteList",object))
})

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #

setGeneric("rewind",function(object)standardGeneric("rewind"))
setMethod("rewind","bamReader",function(object) {return(invisible(.Call("bam_reader_seek",object@reader,object@startpos,PACKAGE="rbamtools")))})

setMethod("bamSave","bamReader",function(object,writer){
  if(!is(writer,"bamWriter"))
    stop("[bamSave.bamReader] writer must be 'bamWriter'!")
  if(!isOpen(object))
    stop("[bamSave.bamReader] reader is not open! Check 'isOpen'!")
  if(!isOpen(writer))
    stop("[bamSave.bamReader] writer is not open! Check 'isOpen'!")
  
  # Saving old reading position
  oldpos<-.Call("bam_reader_tell",object@reader,PACKAGE="rbamtools")
  # Reset reader to start position
  .Call("bam_reader_seek",object@reader,object@startpos,PACKAGE="rbamtools")
  
  nAligns<-.Call("bam_reader_save_aligns",object@reader,writer@writer,PACKAGE="rbamtools")
  bm<-Sys.localeconv()[7]
  cat("[bamSave.bamReader] Saving ",format(nAligns,big.mark=bm)," to file '",basename(writer@filename),"' finished.\n",sep="")
  .Call("bam_reader_seek",object@reader,oldpos,PACKAGE="rbamtools")
  return(invisible(nAligns))
})

setGeneric("bamCopy",function(object,writer,refids,verbose=FALSE)standardGeneric("bamCopy"))
setMethod("bamCopy","bamReader",function(object,writer,refids,verbose=FALSE)
{
  if(!is(writer,"bamWriter"))
    stop("[bamCopy.bamReader] writer must be 'bamWriter'!")
  if(!isOpen(object))
    stop("[bamCopy.bamReader] reader is not open! Check 'isOpen'!")
  if(!isOpen(writer))
    stop("[bamCopy.bamReader] writer is not open! Check 'isOpen'!")
  if(!index.initialized(object))
    stop("[bamCopy.bamReader] reader must have initialized index! Check 'index.initialized'!")
  
  # Check refids argument: When missing copy all ref's
  ref<-getRefData(object)
  if(missing(refids))
  {
    refids<-ref$ID
    n<-length(refids)
    mtc<-1:n
  }
  else
  {
    mtc<-match(refids,ref$ID)
    if(any(is.na(mtc)))
      stop("[bamCopy.bamReader] refids must be subset of Reference-ID's! Check 'getRefData'!")
    n<-length(refids)    
  }

  # Copy aligns with bamRanges as intermediate buffer
  bm<-Sys.localeconv()[7]
  nAligns<-0
  for(i in 1:n)
  {
    range<-bamRange(object,c(ref$ID[mtc[i]],0,ref$LN[mtc[i]]),complex=FALSE)
    nAligns<-nAligns+size(range)
    if(verbose)
      cat("[bamCopy.bamReader] i: ",i,"\tCopying ",format(size(range),big.mark=bm,width=10)," aligns for Reference '",ref$SN[mtc[i]],"'.\n",sep="")
    bamSave(writer,range)
    rm(range)
    gc()    
  }
  cat("[bamCopy.bamReader] Copying ",format(nAligns,big.mark=bm,width=10)," aligns finished.\n",sep="")
})

###################################################################################################
#                                                                                                 #
# bamHeader                                                                                       #
# Description: See SAM File Format Specification (v1.4-r985) September 7,2011, Section 1.3        #
#                                                                                                 #
###################################################################################################

setClass("bamHeader",representation(header="externalptr"),
         validity=function(object){return(ifelse(is.null(object@header),FALSE,TRUE))})

setMethod("initialize","bamHeader",function(.Object,extptr){
  if(!is(extptr,"externalptr"))
    stop("[initialize.bamHeader] extptr must be externalptr!")
  .Object@header<-extptr
  return(.Object)
})

setMethod(f="getHeaderText",signature="bamHeader",definition=function(object) {
  return(new("bamHeaderText",.Call("bam_header_get_header_text",object@header,PACKAGE="rbamtools"))) })

setMethod("as.character","bamHeader",function(x,...){
  .Call("bam_header_get_header_text",x@header,PACKAGE="rbamtools")
})

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #
#  This is the main function for creating an instance of bamWriter              #

setGeneric("bamWriter",function(x,filename)standardGeneric("bamWriter"))
setMethod("bamWriter","bamHeader",function(x,filename){
  if(!is.character(filename))
    stop("[bamWriter.bamHeader] filename must be character!")
  return(new("bamWriter",x,filename))
})
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #


# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #
# headerLine: Represents two entries: Format version (VN) and sorting order (SO)
# Valid format for VN : /^[0-9]+\.[0-9]+$/.
# Valid entries for SO: unknown (default), unsorted, queryname, coordinate.

setClass("headerLine",representation(VN="character",SO="character"),
         validity=function(object)
         {
           if(length(VN)==1 & length(SO)==1)
             return(TRUE)
           else
             return(FALSE)
         })

setMethod(f="initialize",signature="headerLine",definition=function(.Object,hl="",delim="\t"){
  # Parses header line from header section
  if(!is.character(hl))
    stop("[headerLine.initialize] Argument must be string.\n")
  
  # Default object content (hl="" or character(0))
  if((length(hl)==1 && nchar(hl)==0) || length(hl)==0)
  {
    .Object@VN="1.4"
    .Object@SO="unknown"
    return(.Object)
  }
  
  # +++++++++++++++++++++++++++++++++++++++++++++++++++++++
  #  Split input string into tags
  tags<-unlist(strsplit(hl,delim))
  #  Three tags!
  if(length(tags)!=3)
    stop("[headerLine.initialize] hl must contain three tags separated by '",delim,"'!\n")
  #  First  tag: '@HD'
  if(tags[1]!="@HD")
    stop("[headerLine.initialize] First tag of string must be @HD!\n")
  #  Second tag: 'VN'
  #  TODO: Check Accepted format: /^[0-9]+\.[0-9]+$/.  
  if(substr(tags[2],1,2)!="VN")
    stop("[headerLine.initialize] Second tag of string must be VN!\n")
  .Object@VN=substring(tags[2],4)
  # Third   tag: 'SO'
  if(substr(tags[3],1,2)!="SO")
    stop("[headerLine.initialize] Third tag of string must be SO!\n")
  
  str<-substring(tags[3],4)
  if(str=="coordinate")
    .Object@SO<-"coordinate"
  else if(str=="unknown")
    .Object@SO<-"unknown"
  else if(str=="unsorted")
    .Object@SO<-"unsorted"
  else if(str=="queryname")
    .Object@SO<-"queryname"
  
  return(.Object)
})

setMethod("getHeaderText","headerLine",function(object,delim="\t")
{return(paste("@HD\tVN:",object@VN,"\tSO:",object@SO,sep=""))})

setMethod("getVal",signature="headerLine",definition=function(object,member){
  if(!is.character(member))
    stop("[getVal.headerLine] Member must be character!\n")
  if(member=="VN")
    return(object@VN)
  if(member=="SO")
    return(object@SO)
  stop("[getVal.headerLine] Member '",member,"' must be 'VN' or 'SO'!\n")
})

setMethod("setVal",signature="headerLine",definition=function(object,members,values){
  if(!is.character(members) || !is.character(values))
    stop("[setVal.headerLine] Members and values must be character!\n")
  if(length(members)!=length(values))
    stop("[setVal.headerLine] Members and values must have same length!\n")
  tagLabs<-c("VN","SO")
  mtc<-match(members,tagLabs)
  if(any(is.na(mtc)))
    stop("[setVal.headerLine] Member names must be valid Header line entries!\n")
  n<-length(members)
  if(n>2)
    stop("[setVal.headerLine] Only two members can be set!\n")
  obj<-deparse(substitute(object))
  for(i in 1:n)
  {
    txt<-paste(obj,"@",members[i],"<-'",values[i],"'",sep="")
    eval.parent(parse(text=txt))
  }
  return(invisible())
})

setMethod("as.list",signature="headerLine",definition=function(x,...)
{return(list(VN=x@VN,SO=x@SO))})

setMethod("show","headerLine",function(object)
{
  cat("An object of class \"",class(object),"\"\n",sep="")
  cat("VN: ",object@VN,"\nSO: ",object@SO,"\n",sep="")
})

#  End headerLine
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #
#  refSeqDict: Reference Sequence Dictionary                                    #
#  Represents a variable number of Ref Seqs                                     #
#  Valid Members (Entries for each sequence, stored in a data.frame):           #
#  SN Reference sequence name                                                   #
#  LN Reference sequence length                                                 #
#  AS Genome assembly identifier                                                #
#  M5 MD5 checksum of the sequence                                              #
#  SP Species                                                                   #
#  UR URI of the sequence                                                       #
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #

setClass("refSeqDict",representation(SN="character",LN="numeric",AS="character",M5="numeric",SP="character",UR="character"))
setMethod(f="initialize",signature="refSeqDict",definition=function(.Object,hsq="",delim="\t")
{
  # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #
  # Parses Reference sequence dictionary of header-text
  # hsq= Vector of characters, each representing one Ref-Sequence
  # length(hsq) = number of Ref-Sequences
  # Each Ref-string contains 'internally' [tab] delimited seqments:
  #               "SN:ab\tLN:12\tAS:ab\tM5:12\tSP:ab\tUR:ab"
  # It's allowed to skip segments
  # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #
  
  if(!is.character(hsq))
    stop("[refSeqDict.initialize] hsq must be character!")
  
  n<-length(hsq)
  # Return empty object when no input string is given
  if((n==1 && nchar(hsq)==0) || n==0)
    return(.Object)
  
  .Object@SN<-character(n)
  .Object@LN<-numeric(n)
  .Object@AS<-character(n)
  .Object@M5<-numeric(n)
  .Object@SP<-character(n)
  .Object@UR<-character(n)
  labels<-c("SN","LN","AS","M5","SP","UR")
  for(i in 1:n)
  {
    # Containes separated tags for one sequence
    seq<-unlist(strsplit(hsq[i],delim))
    if(seq[1]!="@SQ")
      stop("[initialize.refSeqDict] First segment in Ref-sequence tag must be '@SQ'!")
    seq<-seq[-1]
    
    # Contains column number in dict@df for each tag
    cols<-match(substr(seq,1,2),labels)
    m<-length(cols)
    for(j in 1:m)
    {
      txt<-substr(seq[j],4,nchar(seq[j]))
      # Empty entries are skipped (to avoid errors)
      if(nchar(txt)>0)
      {
        if(cols[j]==1)
          .Object@SN[i]<-txt
        else if(cols[j]==2)
        {
          # Try to convert into numeric value
          numb<-suppressWarnings(as.numeric(txt))  
          if(is.na(numb))
          {warning("[refSeqDict.initialize] No numeric value for LN: '",txt,"'!\n",sep="")}
          else
            .Object@LN[i]<-numb
        }
        else if(cols[j]==3)
          .Object@AS[i]<-txt
        else if(cols[j]==4)
        {
          # Try to convert into numeric value
          numb<-suppressWarnings(as.numeric(txt))  
          if(is.na(numb))
          {warning("[refSeqDict.initialize] No numeric value for LN: '",txt,"'!\n",sep="")}
          else
            .Object@M5<-numb
        }
        else if(cols[j]==5)
          .Object@SP[i]<-txt
        else if(cols[j]==6)
          .Object@UR[i]<-txt
      }
    }
  }
  return(.Object)
})

setMethod(f= "[",signature="refSeqDict",definition=function(x,i){
  rsd<-new("refSeqDict")
  rsd@SN<-x@SN[i]
  rsd@LN<-x@LN[i]
  rsd@AS<-x@AS[i]
  rsd@M5<-x@M5[i]
  rsd@SP<-x@SP[i]
  rsd@UR<-x@UR[i]
  return(rsd)
})

setMethod(f="dim",signature="refSeqDict",definition=function(x){return(c(length(x@SN),6))})

setGeneric("removeSeqs",function(x,rows)standardGeneric("removeSeqs"))
setMethod("removeSeqs",signature="refSeqDict",definition=function(x,rows){
  # Removes given rows (=Sequences) from Dictionary so they are excluded from header
  n<-length(x@SN)
  if(!is.numeric(rows))  
    stop("[removeSeqs.refSeqDict] Sequence indices must be numeric!")
  rows<-as.integer(rows)
  if(any(rows)<1)
    stop("[removeSeqs.refSeqDict] Sequence indices must be positive!")
  if(any(rows)>n)
    stop("[removeSeqs.refSeqDict] Sequence indices must be <",n,"!")
  
  # Execute per eval in parent.frame
  if(length(rows)>1)
    rmv<-paste("c(",paste(rows,collapse=","),")",sep="")
  else
    rmv<-rows
  obj<-deparse(substitute(x))
  dictcol<-paste(obj,"@SN",sep="")
  eval.parent(parse(text=paste(dictcol,"<-",dictcol,"[-",rmv,"]",sep="")))
  
  dictcol<-paste(obj,"@LN",sep="")
  eval.parent(parse(text=paste(dictcol,"<-",dictcol,"[-",rmv,"]",sep="")))
  
  dictcol<-paste(obj,"@AS",sep="")
  eval.parent(parse(text=paste(dictcol,"<-",dictcol,"[-",rmv,"]",sep="")))
  
  dictcol<-paste(obj,"@M5",sep="")
  eval.parent(parse(text=paste(dictcol,"<-",dictcol,"[-",rmv,"]",sep="")))
  
  dictcol<-paste(obj,"@SP",sep="")
  eval.parent(parse(text=paste(dictcol,"<-",dictcol,"[-",rmv,"]",sep="")))
  
  dictcol<-paste(obj,"@UR",sep="")
  eval.parent(parse(text=paste(dictcol,"<-",dictcol,"[-",rmv,"]",sep="")))
  return(invisible())
})

setGeneric("addSeq",function(object,SN,LN,AS="",M5=0,SP="",UR="")standardGeneric("addSeq"))
setMethod("addSeq",signature="refSeqDict",definition=function(object,SN,LN,AS="",M5=0,SP="",UR=""){
  index<-length(object@SN)+1
  obj<-deparse(substitute(object))
  colidx<-paste("[",index,"]",sep="")
  
  # Appends new Sequence (row) at the end
  dictcol<-paste(obj,"@SN",colidx,sep="")
  eval.parent(parse(text=paste(dictcol,"<-'",SN,"'",sep="")))
  
  dictcol<-paste(obj,"@LN",colidx,sep="") 
  eval.parent(parse(text=paste(dictcol,"<-",LN,sep="")))
  
  dictcol<-paste(obj,"@AS",colidx,sep="") 
  eval.parent(parse(text=paste(dictcol,"<-'",AS,"'",sep="")))
  
  dictcol<-paste(obj,"@M5",colidx,sep="") 
  eval.parent(parse(text=paste(dictcol,"<-",M5,sep="")))
  
  dictcol<-paste(obj,"@SP",colidx,sep="") 
  eval.parent(parse(text=paste(dictcol,"<-'",SP,"'",sep="")))
  
  dictcol<-paste(obj,"@UR",colidx,sep="") 
  eval.parent(parse(text=paste(dictcol,"<-'",UR,"'",sep="")))
  
  return(invisible())
})

setMethod("getHeaderText",signature="refSeqDict",definition=function(object,delim="\t"){ 
  # Returns Ref Data String (can be used for creating new BAM file via bamWriter)
  labels<-c("SN","LN","AS","M5","SP","UR")
  n<-length(object@SN)
  if(n==0)
    return(character(0))
  
  seqs<-character(n)
  
  for(i in 1:n)
  {
    ans<-"@SQ"    
    if(nchar(object@SN[i])>0)
      ans<-paste(ans,delim,"SN:",object@SN[i],sep="")
    if(object@LN[i]>0)
      ans<-paste(ans,delim,"LN:",object@LN[i],sep="")
    if(nchar(object@AS[i])>0)
      ans<-paste(ans,delim,"AS:",object@AS[i],sep="")
    if(object@M5[i]>0)
      ans<-paste(ans,delim,"M5:",object@M5[i],sep="")
    if(nchar(object@SP[i])>0)
      ans<-paste(ans,delim,"SP:",object@SP[i],sep="")
    if(nchar(object@UR[i])>0)
      ans<-paste(ans,delim,"UR:",object@UR[i],sep="")
    seqs[i]<-ans
  }
  return(paste(seqs,collapse="\n"))
})

# Return first or last part of refSeqDict data.frame
# S3 Generic is supplied via importFrom in NAMESPACE
setGeneric("head",function(x,...) standardGeneric("head"))
setMethod("head","refSeqDict",function(x,n=6L,...) {
  stopifnot(length(n) == 1L)
  if (n < 0L)
    stop("[head.refSeqDict] n<0!")
  m<-length(x@SN)
  if(m==0)
    cat("[head.refSeqDict] Empty object.\n")
  
  n<-min(n,m)
  if(n == 0L)
    return(as.data.frame(new("refSeqDict")))
  else
    return(as.data.frame(x)[1:n,])
})
# S3 Generic is supplied via importFrom in NAMESPACE
setGeneric("tail",function(x,...) standardGeneric("tail"))
setMethod("tail","refSeqDict",definition=function(x,n=6L,...) {
  stopifnot(length(n) == 1L)
  if (n < 0L)
    stop("[tail.refSeqDict] n<0!")
  m<-length(x@SN)
  if(m==0)
    cat("[tail.refSeqDict] Empty object.\n") 
  n<-min(n,m)
  if(n == 0L)
    return(as.data.frame(new("refSeqDict")))
  else
  {
    n<-m-n+1
    return(x@df[n:m,])
  }
})

setMethod("show","refSeqDict",function(object){
  cat("An object of class \"",class(object),"\"\n",sep="")
  print(head(object))
})

#  End refSeqDict
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #


# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #
# headerReadGroup
# ReadGroup
# ID Read Group identifier
# CN Name of sequencing center
# DS Description
# FO Flow order
# KS Nucleotides corresponding to key sequence of each read
# LB Library
# PG Programs used for processing the Read Group
# PI Predicted median insert size
# PL Sequencing Platform:
#    CAPILLARY,LS454,ILLUMINA,SOLID,HELICOS,IONTORRENT or PACBIO
# SM Sample name.


setClass("headerReadGroup",representation(l="list"),validity=function(object) {return(TRUE)})

setMethod(f="initialize",signature="headerReadGroup", definition=function(.Object,hrg="",delim="\t"){
  # Parses Read-Group part of Header data. See Sam Format Specificatioin 1.3 (Header Section)
  .Object@l<-list()
  if(!is.character(hrg))
    stop("[headerReadGroup.initialize] Argument must be string.\n")
  # hrg="" or character(0)
  if((length(hrg)==1 && nchar(hrg)==0) || length(hrg)==0)
    return(.Object)
  # Split string into fields
  tags<-unlist(strsplit(hrg,delim))
  if(tags[1]!="@RG")
    stop("[headerReadGroup.initialize] First item of string must be @RG!\n")
  tags<-tags[-1]
  tagLabs<-c("ID","CN","DS","DT","FO","KS","LB","PG","PI","PL","PU","SM")
  n<-length(tags)
  for(i in 1:n)
  {
    f<-substr(tags[i],1,2)
    mtc<-match(f,tagLabs)
    if(is.na(mtc))
      stop("[headerReadGroup.initialize] Field identifier '",f,"' not in List!\n")
    .Object@l[[f]]<-substr(tags[i],4,nchar(tags[i]))
  }
  return(.Object)
})

setMethod("getHeaderText",signature="headerReadGroup",definition=function(object,delim="\t") {
  n<-length(object@l)
  if(n==0)
    return(character(0))
  rfstr<-character(n)
  for(i in 1:n)
    rfstr[i]<-paste(names(object@l)[i],object@l[[i]],sep=":")
  return(paste("@RG",paste(rfstr,collapse=delim),sep=delim))
})

setMethod("getVal",signature="headerReadGroup",definition=function(object,member) {
  if(!is.character(member))
    stop("[getVal.headerReadGroup] Member must be character!\n")
  tagLabs<-c("ID","CN","DS","DT","FO","KS","LB","PG","PI","PL","PU","SM")
  mtc<-match(member[1],tagLabs)
  if(is.na(mtc))
    stop("[getVal.headerReadGroup] Invalid member name!\n")
  return(object@l[[member]])
})

setMethod("setVal",signature="headerReadGroup",definition=function(object,members,values){
  if(!is.character(members) || !is.character(values))
    stop("[setVal.headerReadGroup] Member name and value must be character!\n")
  if(length(members)!=length(values))
    stop("[setVal.headerReadGroup] members and values must have same length!\n")
  tagLabs<-c("ID","CN","DS","DT","FO","KS","LB","PG","PI","PL","PU","SM")
  mtc<-match(members,tagLabs)
  if(any(is.na(mtc)))
    stop("[setVal.headerReadGroup] Members must be valid Read Group Entries (See SAM Format Specification 1.3!\n")
  n<-length(members)
  obj<-deparse(substitute(object))
  for(i in 1:n)
  {
    txt<-paste(obj,"@l$",members[i],"<-'",values[i],"'",sep="")
    eval.parent(parse(text=txt))
  }
  return(invisible())
})

setMethod("as.list",signature="headerReadGroup",definition=function(x,...){return(x@l)})

#  End headerReadGroup
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #


# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #
# headerProgram
setClass("headerProgram",representation(l="list"),validity=function(object){return(TRUE)})

setMethod(f="initialize",signature="headerProgram",
          definition=function(.Object,hp="",delim="\t")
          {
            # Parses Program part of Header data. See Sam Format Specificatioin 1.3 (Header Section)
            .Object@l<-list()
            if(!is.character(hp))
              stop("[headerProgram.initialize] Argument must be string.\n")
            # hp="" or character(0)
            if((length(hp)==1 && nchar(hp)==0)||length(hp)==0)
              return(.Object)
            # Split string into fields
            tags<-unlist(strsplit(hp,delim))
            if(tags[1]!="@PG")
              stop("[headerProgram.initialize] First item of string must be @PG!\n")
            tags<-tags[-1]
            tagLabs<-c("ID","PN","CL","PP","VN")
            n<-length(tags)
            for(i in 1:n)
            {
              f<-substr(tags[i],1,2)
              mtc<-match(f,tagLabs)
              if(is.na(mtc))
                stop("[heaProgram.initialize] Field identifier '",f,"' not in List!\n")
              .Object@l[[f]]<-substr(tags[i],4,nchar(tags[i]))
            }
            return(.Object)
          })

setMethod("getHeaderText",signature="headerProgram",definition=function(object,delim="\t") {
  n<-length(object@l)
  if(n==0)
    return(character(0))
  
  rfstr<-character(n)
  for(i in 1:n)
    rfstr[i]<-paste(names(object@l)[i],object@l[[i]],sep=":")
  return(paste("@PG",paste(rfstr,collapse=delim),sep=delim))
})

setMethod("getVal",signature="headerProgram",definition=function(object,member) {
  if(!is.character(member))
    stop("[getVal.headerProgram] Member must be character!\n")
  tagLabs<-c("ID","PN","CL","PP","VN")
  mtc<-match(member[1],tagLabs)
  if(is.na(mtc))
    stop("[getVal.headerProgram] Invalid member name!\n")
  return(object@l[[member]])
})

setMethod("setVal",signature="headerProgram",definition=function(object,members,values){
  if(!is.character(members) || !is.character(values))
    stop("[setVal.headerProgram] Member name and value must be character!\n")
  if(length(members)!=length(values))
    stop("[setVal.headerProgram] members and values must have same length!\n")
  tagLabs<-c("ID","PN","CL","PP","VN")
  mtc<-match(members,tagLabs)
  if(any(is.na(mtc)))
    stop("[setVal.headerProgram] Members must be valid Program Entries (See SAM Format Specification 1.3!\n")
  n<-length(members)
  obj<-deparse(substitute(object))
  for(i in 1:n)
  {
    txt<-paste(obj,"@l$",members[i],"<-'",values[i],"'",sep="")
    eval.parent(parse(text=txt))
  }
  return(invisible())
})

setMethod("as.list",signature="headerProgram",definition=function(x,...){return(x@l)})

setMethod("show","headerProgram",function(object)
{
  cat("An object of class \"",class(object),"\"\n")
  for(i in 1:length(object@l))
  {
    cat(names(object@l)[i],":",object@l[[i]],"\n")
  }
  return(invisible())
})


#  End headerProgram
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #
#  bamHeaderText: Represents and manages textual version of bamHeader           #
#  See SAM Format Specification (v1.4-r985)                                     #
#                                                                               #
#  Contains header Segments :                                                   #
#   head  = headerLine        : @HD Header Line                                 #
#   dict  = refSeqDict        : @SQ Reference Sequence dictionary               #
#   group = headerReadGroup   : @RG Read Group                                  #
#   prog  = headerProgram     : @PG Program                                     #
#                                                                               #
#   TODO:
#   com   = headerComment     : @CO One-line text comment                       #
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #
#  Class definition and creational routines for bamHeaderText

setClass("bamHeaderText",representation(head="headerLine",dict="refSeqDict",
                                        group="headerReadGroup",prog="headerProgram",com="character"))

setMethod(f="initialize",signature="bamHeaderText", definition=function(.Object,bh="",delim="\n")
{
  # Parses Header data (as reported by getHeaderText)
  # See Sam Format Specification 1.3 (Header Section)
  if(!is.character(bh))
    stop("[bamHeaderText.initialize] Argument must be string.\n")
  
  # Create empty header Set (so it's legal to call getHeaderText()')
  if(length(bh)==1 && nchar(bh)==0)
  {
    .Object@head<-new("headerLine")
    .Object@dict<-new("refSeqDict")
    .Object@group<-new("headerReadGroup")
    .Object@prog<-new("headerProgram")
    return(.Object)
  }
  
  # Split input string: Each fragment contains data for one header segment
  bht<-unlist(strsplit(bh,split=delim))
  
  # Read Header Line
  bhl<-bht[grep("@HD",bht)]
  .Object@head<-new("headerLine",bhl)
  
  # Read Sequence Directory
  bsd<-bht[grep("@SQ",bht)]
  .Object@dict<-new("refSeqDict",bsd)
  
  # Read Group
  brg<-bht[grep("@RG",bht)]
  .Object@group<-new("headerReadGroup",brg)
  
  # Read Program Data
  bpd<-bht[grep("@PG",bht)]
  .Object@prog<-new("headerProgram",bpd)
  
  # Read Text comment
  btc<-bht[grep("@CO",bht)]
  com<-substring(btc,3)
  return(.Object)
})

bamHeaderText<-function(head=NULL,dict=NULL,group=NULL,prog=NULL,com=NULL)
{
  bh<-new("bamHeaderText")
  if(!is.null(head))
  {
    if(is(head,"headerLine"))
      bh@head<-head
    else
      stop("[bamHeaderText] head must be 'headerLine'!")
  }
  if(!is.null(dict))
  {
    if(is(dict,"refSeqDict"))
      bh@dict<-dict
    else
      stop("[bamHeaderText] dict must be 'refSeqDict'")
  }
  if(!is.null(group))
  {
    if(is(group,"headerReadGroup"))
      bh@group<-group
    else
      stop("[bamHeaderText] group must be 'headerReadGroup'!")
  }
  if(!is.null(prog))
  {
    if(is(prog,"headerProgram"))
      bh@prog<-prog
    else
      stop("[bamHeaderText] prog must be 'headerProgram'!")
  }
  if(!is.null(com))
  {
    if(is.character(com))
      bh@com<-com
    else
      stop("[bamHeaderText] com must be 'character'!")
  }
  return(invisible(bh))
}

#  End: Class definition and creational routines for bamHeaderText
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #
#  Public accessors for member objects for bamHeaderText
setGeneric("headerLine",function(object) standardGeneric("headerLine"))
setMethod(f="headerLine",signature="bamHeaderText",definition=function(object) {return(object@head)})
setGeneric("refSeqDict",function(object) standardGeneric("refSeqDict"))
setMethod(f="refSeqDict",signature="bamHeaderText",definition=function(object) {return(object@dict)})
setGeneric("headerReadGroup",function(object)standardGeneric("headerReadGroup"))
setMethod(f="headerReadGroup",signature="bamHeaderText",definition=function(object){return(object@group)})
setGeneric("headerProgram",function(object)standardGeneric("headerProgram"))
setMethod(f="headerProgram",signature="bamHeaderText",definition=function(object){return(object@prog)})

setGeneric("headerLine<-",function(object,value)standardGeneric("headerLine<-"))
setReplaceMethod("headerLine","bamHeaderText",function(object,value)
{
  if(!is(value,"headerLine"))
    stop("[headerLine<-.bamHeaderText] value must be 'headerLine'!")
  object@head<-value
  return(object)
})

setGeneric("refSeqDict<-",function(object,value)standardGeneric("refSeqDict<-"))
setReplaceMethod("refSeqDict","bamHeaderText",function(object,value)
{
  if(!is(value,"refSeqDict"))
    stop("[refSeqDict<-.bamHeaderText] value must be 'refSeqDict'!")
  object@dict<-value
  return(object)
})

setGeneric("headerReadGroup<-",function(object,value)standardGeneric("headerReadGroup<-"))
setReplaceMethod("headerReadGroup","bamHeaderText",function(object,value)
{
  if(!is(value,"headerReadGroup"))
    stop("[headerReadGroup<-.bamHeaderText] value must be 'headerReadGroup'!")
  object@group<-value
  return(object)
})

setGeneric("headerProgram<-",function(object,value)standardGeneric("headerProgram<-"))
setReplaceMethod("headerProgram","bamHeaderText",function(object,value)
{
  if(!is(value,"headerProgram"))
    stop("[headerProgram<-.bamHeaderText] value must be 'headerProgram'!")
  object@prog<-value
  return(object)
})

#  End: Public accessors for member objects for bamHeaderText
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #



setMethod("getHeaderText",signature="bamHeaderText",definition=function(object,delim="\n") {
  hd<-getHeaderText(object@head)
  if(length(hd)==0)
    return(character(0))
  hd<-paste(hd,delim,sep="")
  
  dt<-getHeaderText(object@dict)
  if(length(dt)==0)
    return(character(0))
  dt<-paste(dt,delim,sep="")
  
  gp<-getHeaderText(object@group)
  if(length(gp)>0)
    gp<-paste(gp,delim,sep="")
  
  pg<-getHeaderText(object@prog)
  if(length(pg)>0)
    pg<-paste(pg,delim,sep="")
  
  if(length(object@com)>0)
    cm<-paste(paste("@CO",object@com,sep="\t"),collapse=delim)
  else
    cm<-character(0)
  return(paste(hd,dt,gp,pg,cm,sep=""))
})

setGeneric("bamHeader",function(object)standardGeneric("bamHeader"))
setMethod("bamHeader","bamHeaderText",
          function(object){return(new("bamHeader",.Call("init_bam_header",getHeaderText(object))))})


###################################################################################################
#                                                                                                 #
# bamWriter class                                                                                 #
# Encapsulates an write-opened Connection to a BAM-file.                                          #
#                                                                                                 #
###################################################################################################

setClass("bamWriter",representation(filename="character",writer="externalptr"),
         validity=function(object) {return(ifelse(is.null(object@writer),FALSE,TRUE))})

setMethod(f="initialize", signature="bamWriter",
          definition=function(.Object,header,filename){
            if(!is(header,"bamHeader"))
              stop("[initialize.bamWriter] header must be bamHeader!\n")
            if(!is.character(filename))
              stop("[initialize.bamWriter] filename must be character!\n")
            .Object@filename<-filename
            .Object@writer<-.Call("bam_writer_open",header@header,filename,PACKAGE="rbamtools")
            return(.Object)
          })


setMethod("filename", "bamWriter", function(object) return(object@filename))
setMethod("isOpen",signature="bamWriter",definition=function(con,rw="")
{return(!(.Call("is_nil_externalptr",con@writer,PACKAGE="rbamtools")))})

setMethod(f="bamClose",signature="bamWriter",definition=function(object)
{ invisible(.Call("bam_writer_close",object@writer,PACKAGE="rbamtools"))})

setMethod(f="bamSave",signature="bamWriter",definition=function(object,value) 
{
  if(is(value,"bamAlign"))
  {
    return(invisible(.Call("bam_writer_save_align",object@writer,value@align,PACKAGE="rbamtools")))
  }
  if(is(value,"bamRange"))
  {
    return(invisible(.Call("bam_range_write",object@writer,value@range,PACKAGE="rbamtools")))
  }
  else
    stop("bamSave: Saved object must be of type bamAlign or bamRange!\n")
})


###################################################################################################
#                                                                                                 #
# gapList                                                                                         #
#                                                                                                 #
###################################################################################################

setClass("gapList",representation(list="externalptr"),
         validity=function(object){ return(ifelse(is.null(object@list),FALSE,TRUE))})

setMethod(f="initialize","gapList",definition=function(.Object,reader,coords,verbose=FALSE){
  
  if(!is(reader,"bamReader"))
  {
    cat("[initialize.gapList] Class of reader: ",class(reader),".\n")
    stop("[initialize.gapList] reader must be an instance of bamReader!\n")
  }
  if(length(coords)!=3)
    stop("[initialize.gapList] coords must be 3-dim numeric (ref,start,stop)!\n")  
  if(is.null(reader@index))
    stop("[initialize.gapList] bamReader must have initialized index!\n")
  .Object@list<-.Call("gap_list_fetch",reader@reader,reader@index,trunc(coords),PACKAGE="rbamtools")
  glsize<-.Call("gap_list_get_size",.Object@list,PACKAGE="rbamtools")
  if(verbose)
    message("[initialize.gapList] Fetched list of size ",format(glsize,big.mark=Sys.localeconv()[7])," for refid ",coords[1],".")
  return(.Object)
})

# gapList function for retrieving objects in bamReader section

setMethod("size",signature="gapList",definition=function(object)
{.Call("gap_list_get_size",object@list,PACKAGE="rbamtools")})

setMethod("nAligns",signature="gapList",definition=function(object)
{.Call("gap_list_get_nAligns",object@list,PACKAGE="rbamtools")})

setMethod("nGapAligns",signature="gapList",definition=function(object)
{.Call("gap_list_get_nGapAligns",object@list,PACKAGE="rbamtools")})

setMethod("show","gapList",function(object){
  cat("An object of class '",class(object),"'. size: ",size(object),"\n",sep="")
  cat("nAligns:",nAligns(object),"\tnGapAligns:",nGapAligns(object),"\n")
  return(invisible())
})

###################################################################################################
#                                                                                                 #
# gapSiteList                                                                                     #
#                                                                                                 #
###################################################################################################

setClass("gapSiteList",representation(list="externalptr"),
         validity=function(object){ return(ifelse(is.null(object@list),FALSE,TRUE))})

setMethod(f="initialize","gapSiteList",definition=function(.Object,reader,coords){
  if(missing(reader) || missing(coords))
    return(.Object)
  
  if(!is(reader,"bamReader"))
    stop("[gapSiteList.initialize] reader must be an instance of bamReader!\n")
  if(length(coords)!=3)
    stop("[gapSiteList.initialize] coords must be 3-dim numeric (ref,start,stop)!\n")  
  if(is.null(reader@index))
    stop("[gapSiteList.initialize] bamReader must have initialized index!\n")
  .Object@list<-.Call("gap_site_list_fetch",reader@reader,reader@index,trunc(coords),PACKAGE="rbamtools")
  return(.Object)
})


# siteList function for retrieving objects in bamReader section

setMethod("size",signature="gapSiteList",definition=function(object)
{.Call("gap_site_list_get_size",object@list,PACKAGE="rbamtools")})

setMethod("nAligns",signature="gapSiteList",definition=function(object)
{.Call("gap_site_list_get_nAligns",object@list,PACKAGE="rbamtools")})

setMethod("nGapAligns",signature="gapSiteList",definition=function(object)
{.Call("gap_site_list_get_nGapAligns",object@list,PACKAGE="rbamtools")})

setMethod("show","gapSiteList",function(object){
  cat("An object of class '",class(object),"'. size: ",size(object),"\n",sep="")
  cat("nAligns:",nAligns(object),"\tnGapAligns:",nGapAligns(object),"\n")
  return(invisible())
})

merge.gapSiteList<-function(x,y,...)
{
  if(!is(y,"gapSiteList"))
    stop("[merge.gapSiteList] y must be of class 'gapSiteList'!")
  res<-new("gapSiteList")
  res@list<-.Call("gap_site_list_merge",x@list,y@list)
  return(res)
}

###################################################################################################
#                                                                                                 #
# bamSiteList                                                                                     #
#                                                                                                 #
###################################################################################################

setClass("bamSiteList",representation(list="externalptr",refdata="data.frame"),
         validity=function(object){ return(ifelse(is.null(object@list),FALSE,TRUE))})

setMethod(f="initialize","bamSiteList",definition=function(.Object,reader){
  if(missing(reader))
  {
    .Object@list<-.Call("gap_site_ll_init")
    return(.Object)
  }
  
  if(!is(reader,"bamReader"))
    stop("[bamSiteList.initialize] reader must be an instance of bamReader!\n")
  if(is.null(reader@index))
    stop("[bamSiteList.initialize] bamReader must have initialized index!\n")
  
  ref<-getRefData(reader)
  ref$start<-0L
  .Object@list<-.Call("gap_site_ll_fetch",reader@reader,reader@index,ref$ID,ref$start,ref$LN,PACKAGE="rbamtools")
  
  # filter refdata for existing lists
  sm<-.Call("gap_site_ll_get_summary_df",.Object@list,PACKAGE="rbamtools")
  mtc<-match(ref$ID,sm$ID)
  .Object@refdata<-ref[!is.na(mtc),]
  
  # Re-enumerate ID's to 1:n
  .Object@refdata$ID<-.Call("gap_site_ll_reset_refid",.Object@list,PACKAGE="rbamtools")
  # ToDo: merge refdata with summary df?
  
  return(.Object)
})


# siteList function for retrieving objects in bamReader section

setMethod("size",signature="bamSiteList",definition=function(object)
{.Call("gap_site_ll_get_size",object@list,PACKAGE="rbamtools")})

setMethod("nAligns",signature="bamSiteList",definition=function(object)
{.Call("gap_site_ll_get_nAligns",object@list,PACKAGE="rbamtools")})

setMethod("nGapAligns",signature="bamSiteList",definition=function(object)
{.Call("gap_site_ll_get_nGapAligns",object@list,PACKAGE="rbamtools")})

setMethod("show","bamSiteList",function(object){
  bm<-Sys.localeconv()[7]
  cat("An object of class '",class(object),"'. size: ",format(size(object),big.mark=bm),"\n",sep="")
  cat("nAligns:",format(nAligns(object),big.mark=bm),"\tnGapAligns:",format(nGapAligns(object),big.mark=bm),"\n")
  
  return(invisible())
})

summary.bamSiteList<-function(object,...)
{ return(merge(object@refdata,.Call("gap_site_ll_get_summary_df",object@list))) }

merge.bamSiteList<-function(x,y,...)
{
  if(!is(y,"bamSiteList"))
    stop("[merge.bamSiteList] y must be bamSiteList!")
  #lref<-x@refdata
  #rref<-y@refdata
  mref<-merge(x@refdata,y@refdata,by="SN",all=T)
  
  n<-dim(mref)[1]
  .Call("gap_site_ll_set_curr_first",x@list)
  .Call("gap_site_ll_set_curr_first",y@list)
  res<-new("bamSiteList")
  for(i in 1:n)
  {
    if(is.na(mref$ID.x[i]))
      .Call("gap_site_ll_add_curr_pp",y@list,res@list,as.integer(i-1))
    else if(is.na(mref$ID.y[i]))
      .Call("gap_site_ll_add_curr_pp",x@list,res@list,as.integer(i-1))
    else
      .Call("gap_site_ll_add_merge_pp",x@list,y@list,res@list,as.integer(i-1))
  }
  
  # get l-part of refdata
  ysn<-mref$SN[is.na(mref$ID.x)]
  mtc<-match(ysn,y@refdata$SN)
  res@refdata<-rbind(x@refdata,y@refdata[mtc,])
  # reset ID to new values
  res@refdata$ID<-0:(n-1)
  return(res)
}

readPooledBamGaps<-function(infiles,idxInfiles=paste(infiles,".bai",sep=""))
{
  bm<-Sys.localeconv()[7]
  n<-length(infiles)
  for(i in 1:n)
  {
    cat("[readPooledBamGaps] (",format(i,width=2),"/",n,") File opened.",sep="")
    bam<-infiles[i]
    reader<-bamReader(bam)
    if(!file.exists(idxInfiles[i]))
      createIndex(reader,idxInfiles[i])
    loadIndex(reader,idxInfiles[i])
    if(i==1)
      ga<-bamSiteList(reader)
    else
    {
      ga1<-bamSiteList(reader)
      ga<-merge.bamSiteList(ga,ga1)
    }
    cat("\r[readPooledBamGaps] (",format(i,width=2),"/",n,") Finished. List-size: ",format(size(ga),width=7,big.mark=bm),
              "\tnAligns: ",format(nAligns(ga),width=13,big.mark=bm),".\n",sep="")
  }
  cat("\n")
  return(ga)
}

readPooledBamGapDf<-function(infiles,idxInfiles=paste(infiles,".bai",sep=""))
{ 
  ga<-readPooledBamGaps(infiles,idxInfiles=paste(infiles,".bai",sep=""))
  dfr<-as.data.frame(ga)
  dfr$gqs<-dfr$nlstart*dfr$qmm
  attr(dfr,"nAligns")<-nAligns(ga)
  attr(dfr,"nGapAligns")<-nGapAligns(ga)
  return(dfr)
}


###################################################################################################
#                                                                                                 #
# bamRange                                                                                        #
#                                                                                                 #
###################################################################################################

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #
# Encapsulates a bunch of Alignment datasets that typically have been read from a defined         #
# reference region in a BAM-file.                                                                 #
# Technically, the alignments are stored in a (C-implemented) double linked list.                 #
# bamRange objects can be created by a reading procedure on an indexed BAM-file. The alignments   #
# can be iterated, readed, written, deleted and added. bamRange objects can be written to a       #
# BAM-file via an Instance of bamWriter.                                                          #
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #

bamRange<-function(reader=NULL,coords=NULL,complex=FALSE) {
  if(!is.null(reader))
  {
    if(!index.initialized(reader))
      stop("[bamRange] reader must have initialized index! Use 'loadIndex'!")
  }
  return(new("bamRange",reader,coords,complex))
}

setClass("bamRange",representation(range="externalptr"),
         validity=function(object) { return(ifelse(is.null(object@range),FALSE,TRUE)) })

setMethod(f="initialize",signature="bamRange",
          definition=function(.Object,reader=NULL,coords=NULL,complex=FALSE)
          {        
            # +++++++++++++++++++++++++++++++++++++++++++            
            #  Create empty range
            if(is.null(reader))
            {
              .Object@range<-.Call("bam_range_init")
              return(.Object)
            }
            
            # +++++++++++++++++++++++++++++++++++++++++++
            #  Create range from bam-file
            if(!is(reader,"bamReader"))
              stop("[bamRange.initialize] reader must be an instance of bamReader!")
            if(length(coords)!=3)
              stop("[bamRange.initialize] coords must be 3-dim numeric (ref,start,stop)!")  
            if(is.null(reader@index))
              stop("[bamRange.initialize] bamReader must have initialized index!")
            if(!is(complex,"logical"))
              stop("[bamRange.initialize] complex must be logical!")
            if(length(complex)>1)
              stop("[bamRange.initialize] complex must have length 1!")
            if(!index.initialized(reader))
              stop("[bamRange.initialize] reader must have initialized index! Use 'loadIndex'!")
            .Object@range<-.Call("bam_range_fetch",reader@reader,reader@index,trunc(coords),complex,PACKAGE="rbamtools")
            return(.Object)
          })

# REMOVED:
# setMethod("as.data.frame",signature="bamRange",definition=function(x,row.names=NULL,optional=FALSE,...) {
#   return(.Call("bam_range_get_align_df",x@range,PACKAGE="rbamtools"))
# })



setMethod("size",signature="bamRange",definition=function(object)
{.Call("bam_range_get_size",object@range,PACKAGE="rbamtools")})

setMethod("show","bamRange",function(object){
  cat("An object of class '",class(object),"'. size: ",size(object),"\n",sep="")
  return(invisible())
})

setMethod("getNextAlign",signature="bamRange",definition=function(object)
{
  ans<-.Call("bam_range_get_next_align",object@range,PACKAGE="rbamtools")
  # Must be checked because align list returns NULL when end is reached
  if(is.null(ans))
    return(ans)
  else
    return(new("bamAlign",ans))
})

setGeneric("getPrevAlign",function(object) standardGeneric("getPrevAlign"))
setMethod("getPrevAlign",signature="bamRange",definition=function(object)
{ return(new("bamAlign",.Call("bam_range_get_prev_align",object@range,PACKAGE="rbamtools")))})

setGeneric("stepNextAlign",function(object)standardGeneric("stepNextAlign"))
setMethod("stepNextAlign",signature("bamRange"),definition=function(object)
{
  .Call("bam_range_step_next_align",object@range)
  return(invisible())
})

setGeneric("stepPrevAlign",function(object)standardGeneric("stepPrevAlign"))
setMethod("stepPrevAlign",signature("bamRange"),definition=function(object)
{
  .Call("bam_range_step_prev_align",object@range)
  return(invisible())
})

# Resets current align to NULL position (i.e. before first element)
# The next call to getNextAlign then returns the first element of list
setGeneric("windBack", function(object) standardGeneric("windBack"))
setMethod("windBack",signature="bamRange",definition=function(object)
{invisible(.Call("bam_range_wind_back",object@range,PACKAGE="rbamtools"))})

setGeneric("push_back",function(object,value) standardGeneric("push_back"))
setMethod("push_back",signature="bamRange",definition=function(object,value)
{
  if(!is(value,"bamAlign"))
    stop("push_back.bamRange: pushed object must be of class \"bamAlign\"\n")
  .Call("bam_range_push_back",object@range,value@align,PACKAGE="rbamtools")
})

setGeneric("pop_back",function(object) standardGeneric("pop_back"))
setMethod("pop_back",signature="bamRange",definition=function(object)
{.Call("bam_range_pop_back",object@range,PACKAGE="rbamtools") })

setGeneric("push_front",function(object,value) standardGeneric("push_front"))
setMethod("push_front",signature="bamRange",definition=function(object,value)
{
  if(!is(value,"bamAlign"))
    stop("push_front.bamRange: pushed object must be of class \"bamAlign\"\n")
  .Call("bam_range_push_front",object@range,value@align,PACKAGE="rbamtools")
})

setGeneric("pop_front",function(object) standardGeneric("pop_front"))
setMethod("pop_front",signature="bamRange",definition=function(object)
{.Call("bam_range_pop_front",object@range,PACKAGE="rbamtools")})

setGeneric("writeCurrentAlign",function(object,value) standardGeneric("writeCurrentAlign"))
setMethod("writeCurrentAlign",signature="bamRange",definition=function(object,value)
{
  if(!is(value,"bamAlign"))
    stop("writeCurrentAlign.bamRange: written object must be of class \"bamAlign\"\n")
  .Call("bam_range_write_current_align",object@range,value@align,PACKAGE="rbamtools")
})

setGeneric("insertPastCurrent",function(object,value) standardGeneric("insertPastCurrent"))
setMethod("insertPastCurrent",signature="bamRange",definition=function(object,value)
{
  if(!is(value,"bamAlign"))
    stop("insertPastCurrent.bamRange: written object must be of class \"bamAlign\"\n")
  .Call("bam_range_insert_past_curr_align",object@range,value@align,PACKAGE="rbamtools")
})

setGeneric("insertPreCurrent",function(object,value) standardGeneric("insertPreCurrent"))
setMethod("insertPreCurrent",signature="bamRange",definition=function(object,value)
{
  if(!is(value,"bamAlign"))
    stop("insertPreCurrent.bamRange: written object must be of class \"bamAlign\"\n")
  .Call("bam_range_insert_pre_curr_align",object@range,value@align,PACKAGE="rbamtools")
})

setGeneric("moveCurrentAlign",function(object,target) standardGeneric("moveCurrentAlign"))
setMethod("moveCurrentAlign",signature="bamRange",definition=function(object,target)
{
  if(!is(target,"bamRange"))
    stop("[moveCurrentAlign.bamRange] target must be bamRange!\n")
  .Call("bam_range_mv_curr_align",object@range,target@range)
  return(invisible())
})

###################################################################################################
#                                                                                                 #
# bamAlign                                                                                        #
#                                                                                                 #
###################################################################################################

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #
# bamAlign encapsulates all contained data in a single dataset in a BAM-file. bamAlign objects    #
# can be read from a bamReader instance and written to a bamWriter instance. All contained data   #
# can be read and written via accessor functions.                                                 #
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #

setClass("bamAlign", representation(align="externalptr"),
         validity=function(object){return(ifelse(is.null(object@align,FALSE,TRUE)))})

setMethod(f="initialize", signature="bamAlign",
          definition=function(.Object,align=NULL){
            .Object@align<-align
            return(.Object)
          }
)


# bamAlign Member Reader functions
setGeneric("name",function(object) standardGeneric("name"))
setMethod(f="name",signature="bamAlign",definition=function(object) 
{ .Call("bam_align_get_name",object@align,PACKAGE="rbamtools") })

setGeneric("refID",function(object) standardGeneric("refID"))
setMethod(f="refID",signature="bamAlign",definition=function(object)
{return(.Call("bam_align_get_refid",object@align,PACKAGE="rbamtools"))})

setGeneric("position",function(object) standardGeneric("position"))
setMethod(f="position",signature="bamAlign",definition=function(object)
{return(.Call("bam_align_get_position",object@align,PACKAGE="rbamtools"))})

setGeneric("nCigar",function(object) standardGeneric("nCigar"))
setMethod("nCigar",signature="bamAlign",definition=function(object)
{ return(.Call("bam_align_get_nCigar",object@align,PACKAGE="rbamtools"))})

setGeneric("cigarData",function(object) standardGeneric("cigarData"))
setMethod(f="cigarData",signature="bamAlign",definition=function(object)
{ .Call("bam_align_get_cigar_df",object@align,PACKAGE="rbamtools")})

setGeneric("mateRefID",function(object) standardGeneric("mateRefID"))
setMethod(f="mateRefID",signature="bamAlign",definition=function(object)
{ .Call("bam_align_get_mate_refid",object@align,PACKAGE="rbamtools")})

setGeneric("matePosition",function(object) standardGeneric("matePosition"))
setMethod(f="matePosition",signature="bamAlign",definition=function(object)
{ .Call("bam_align_get_mate_position",object@align,PACKAGE="rbamtools")})

setGeneric("insertSize",function(object) standardGeneric("insertSize"))
setMethod(f="insertSize",signature="bamAlign",definition=function(object)
{ .Call("bam_align_get_insert_size",object@align,PACKAGE="rbamtools")})

setGeneric("mapQuality",function(object) standardGeneric("mapQuality"))
setMethod(f="mapQuality",signature="bamAlign",definition=function(object)
{ .Call("bam_align_get_map_quality",object@align,PACKAGE="rbamtools")})

setGeneric("alignSeq",function(object) standardGeneric("alignSeq"))
setMethod(f="alignSeq",signature="bamAlign",definition=function(object)
{return(.Call("bam_align_get_segment_sequence",object@align,PACKAGE="rbamtools"))})

setGeneric("alignQual",function(object) standardGeneric("alignQual"))
setMethod(f="alignQual",signature="bamAlign",definition=function(object)
{return(.Call("bam_align_get_qualities",object@align,PACKAGE="rbamtools"))})

setMethod("show","bamAlign",function(object){
  cat("An object of class '",class(object),"'.\n",sep="")
  cat("refID:",refID(object),"\tposition:",position(object),"\n")
  cat("cigarData:\n")
  print(cigarData(object))
})

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #
# Queries against alignment flag (Readers and Accessors)

# pcrORopt_duplicate
setGeneric("pcrORopt_duplicate", function(object) standardGeneric("pcrORopt_duplicate"))
setMethod("pcrORopt_duplicate", "bamAlign", function(object)
  return(.Call("bam_align_is_pcr_or_optical_dup",object@align,PACKAGE="rbamtools")))
setGeneric("pcrORopt_duplicate<-", function(object,value) standardGeneric("pcrORopt_duplicate<-"))
setReplaceMethod(f="pcrORopt_duplicate", signature="bamAlign",
                 definition=function(object,value){
                   if(!is.logical(value))
                     stop("class bamReader, Duplicate setter: value must be boolean")
                   .Call("bam_align_set_is_pcr_or_optical_dup",object@align,value,PACKAGE="rbamtools")
                   return(object)
                 }
)
# failedQC
setGeneric("failedQC", function(object) standardGeneric("failedQC"))
setMethod("failedQC", "bamAlign", function(object)
  return(.Call("bam_align_fail_qc",object@align,PACKAGE="rbamtools")))
setGeneric("failedQC<-", function(object,value) standardGeneric("failedQC<-"))
setReplaceMethod(f="failedQC", signature="bamAlign",
                 definition=function(object,value){
                   if(!is.logical(value))
                     stop("class bamReader, failedQC setter: value must be boolean")
                   .Call("bam_align_set_fail_qc",object@align,value,PACKAGE="rbamtools")
                   return(object)
                 }
)
# firstInPair
setGeneric("firstInPair", function(object) standardGeneric("firstInPair"))
setMethod("firstInPair", "bamAlign", function(object)
  return(.Call("bam_align_is_first_in_pair",object@align,PACKAGE="rbamtools")))
setGeneric("firstInPair<-", function(object,value) standardGeneric("firstInPair<-"))
setReplaceMethod(f="firstInPair", signature="bamAlign",
                 definition=function(object,value){
                   if(!is.logical(value))
                     stop("class bamReader, FirstInPair setter: value must be boolean")
                   .Call("bam_align_set_is_first_in_pair",object@align,value,PACKAGE="rbamtools")
                   return(object)
                 }
)
# secondInPair
setGeneric("secondInPair", function(object) standardGeneric("secondInPair"))
setMethod("secondInPair", "bamAlign", function(object)
  return(.Call("bam_align_is_second_in_pair",object@align,PACKAGE="rbamtools")))
setGeneric("secondInPair<-", function(object,value) standardGeneric("secondInPair<-"))
setReplaceMethod(f="secondInPair", signature="bamAlign",
                 definition=function(object,value){
                   if(!is.logical(value))
                     stop("class bamReader, secondInPair setter: value must be boolean")
                   .Call("bam_align_set_is_second_in_pair",object@align,value,PACKAGE="rbamtools")
                   return(object)
                 }
)
# unmapped
setGeneric("unmapped", function(object) standardGeneric("unmapped"))
setMethod("unmapped", "bamAlign", function(object)
  return(.Call("bam_align_is_unmapped",object@align,PACKAGE="rbamtools")))
setGeneric("unmapped<-", function(object,value) standardGeneric("unmapped<-"))
setReplaceMethod(f="unmapped", signature="bamAlign",
                 definition=function(object,value){
                   if(!is.logical(value))
                     stop("class bamReader, unmapped setter: value must be boolean")
                   .Call("bam_align_set_is_unmapped",object@align,value,PACKAGE="rbamtools")
                   return(object)
                 }
)
# mateUnmapped
setGeneric("mateUnmapped", function(object) standardGeneric("mateUnmapped"))
setMethod("mateUnmapped", "bamAlign", function(object)
  return(.Call("bam_align_mate_is_unmapped",object@align,PACKAGE="rbamtools")))
setGeneric("mateUnmapped<-", function(object,value) standardGeneric("mateUnmapped<-"))
setReplaceMethod(f="mateUnmapped", signature="bamAlign",
                 definition=function(object,value){
                   if(!is.logical(value))
                     stop("class bamReader, mateUnmapped setter: value must be boolean")
                   .Call("bam_align_set_mate_is_unmapped",object@align,value,PACKAGE="rbamtools")
                   return(object)
                 }
)
# reverseStrand
setGeneric("reverseStrand", function(object) standardGeneric("reverseStrand"))
setMethod("reverseStrand", "bamAlign", function(object)
  return(.Call("bam_align_strand_reverse",object@align,PACKAGE="rbamtools")))
setGeneric("reverseStrand<-", function(object,value) standardGeneric("reverseStrand<-"))
setReplaceMethod(f="reverseStrand", signature="bamAlign",
                 definition=function(object,value){
                   if(!is.logical(value))
                     stop("class bamReader, reverseStrand setter: value must be boolean")
                   .Call("bam_align_set_strand_reverse",object@align,value,PACKAGE="rbamtools")
                   return(object)
                 }
)
# mateReverseStrand
setGeneric("mateReverseStrand", function(object) standardGeneric("mateReverseStrand"))
setMethod("mateReverseStrand", "bamAlign", function(object)
  return(.Call("bam_align_mate_strand_reverse",object@align,PACKAGE="rbamtools")))
setGeneric("mateReverseStrand<-", function(object,value) standardGeneric("mateReverseStrand<-"))
setReplaceMethod(f="mateReverseStrand", signature="bamAlign",
                 definition=function(object,value){
                   if(!is.logical(value))
                     stop("class bamReader, mateReverseStrand setter: value must be boolean")
                   .Call("bam_align_set_mate_strand_reverse",object@align,value,PACKAGE="rbamtools")
                   return(object)
                 }
)
# paired
setGeneric("paired", function(object) standardGeneric("paired"))
setMethod("paired", "bamAlign", function(object)
  return(.Call("bam_align_is_paired",object@align,PACKAGE="rbamtools")))
setGeneric("paired<-", function(object,value) standardGeneric("paired<-"))
setReplaceMethod(f="paired", signature="bamAlign",
                 definition=function(object,value){
                   if(!is.logical(value))
                     stop("class bamReader, paired setter: value must be boolean")
                   .Call("bam_align_set_is_paired",object@align,value,PACKAGE="rbamtools")
                   return(object)
                 }
)
# properPair
setGeneric("properPair", function(object) standardGeneric("properPair"))
setMethod("properPair", "bamAlign", function(object)
  return(.Call("bam_align_mapped_in_proper_pair",object@align,PACKAGE="rbamtools")))
setGeneric("properPair<-", function(object,value) standardGeneric("properPair<-"))
setReplaceMethod(f="properPair", signature="bamAlign",
                 definition=function(object,value){
                   if(!is.logical(value))
                     stop("class bamReader, properPair setter: value must be boolean")
                   .Call("bam_align_set_mapped_in_proper_pair",object@align,value,PACKAGE="rbamtools")
                   return(object)
                 }
)
# secondaryAlign
setGeneric("secondaryAlign", function(object) standardGeneric("secondaryAlign"))
setMethod("secondaryAlign", "bamAlign", function(object)
  return(.Call("bam_align_is_secondary_align",object@align,PACKAGE="rbamtools")))
setGeneric("secondaryAlign<-", function(object,value) standardGeneric("secondaryAlign<-"))
setReplaceMethod(f="secondaryAlign", signature="bamAlign",
                 definition=function(object,value){
                   if(!is.logical(value))
                     stop("class bamReader, SecondaryAlign setter: value must be boolean")
                   .Call("bam_align_set_is_secondary_align",object@align,value,PACKAGE="rbamtools")
                   return(object)
                 }
)
# flag
setGeneric("flag", function(object) standardGeneric("flag"))
setMethod("flag", "bamAlign", function(object)
  return(.Call("bam_align_get_flag",object@align,PACKAGE="rbamtools")))
setGeneric("flag<-", function(object,value) standardGeneric("flag<-"))
setReplaceMethod(f="flag", signature="bamAlign",
                 definition=function(object,value){
                   if(!is.integer(value))
                     stop("class bamReader, flag setter: value must be boolean")
                   .Call("bam_align_set_flag",object@align,value,PACKAGE="rbamtools")
                   return(object)
                 }
)
#  End: Queries against alignment flag (Readers and Accessors)
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #

#  End: bamAlign
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #


# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #
#  coercing

as.data.frame.bamRange<-function(x,row.names=NULL,optional=FALSE,...)
  {return(.Call("bam_range_get_align_df",x@range,PACKAGE="rbamtools"))}
as.data.frame.gapList<-function(x,row.names=NULL,optional=FALSE,...)
  {return(.Call("gap_list_get_df",x@list,PACKAGE="rbamtools"))}
as.data.frame.gapSiteList<-function(x,row.names=NULL,optional=FALSE,...)
{return(.Call("gap_site_list_get_df",x@list,PACKAGE="rbamtools"))}
as.data.frame.bamSiteList<-function(x,row.names=NULL,optional=FALSE,...)
{return(.Call("gap_site_ll_get_df",x@list,x@refdata$SN,PACKAGE="rbamtools"))}

as.data.frame.refSeqDict<-function(x,row.names=NULL,optional=FALSE,...)
{
  if(is.null(row.names))
    row.names<-1:(length(x@SN))
  else if(length(row.names)!=length(x@SN))
    stop("[as.data.frame.refSeqDict] length(row.names)!=length(x@SN)!")
  return(data.frame(SN=x@SN,LN=x@LN,AS=x@AS,M5=x@M5,SP=x@SP,UR=x@UR,row.names=row.names))  
}


setAs("bamRange","data.frame",function(from)
  {return(.Call("bam_range_get_align_df",from@range,PACKAGE="rbamtools"))})
setAs("gapList","data.frame",function(from)
  {return(.Call("gap_list_get_df",from@list,PACKAGE="rbamtools"))})
setAs("refSeqDict","data.frame",function(from)
  {return(data.frame(SN=from@SN,LN=from@LN,AS=from@AS,M5=from@M5,SP=from@SP,UR=from@UR,row.names=1:length(from@SN)))})
