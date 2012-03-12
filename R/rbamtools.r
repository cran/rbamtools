
# .First.lib(libname, pkgname) is called when a package is loaded; .Last.lib(libpath) is called when a package is detached.
# exportPattern("^[[:alpha:]]+")
#.First.lib<-function(libname,pkgname) { library.dynam("rbamtools", pkgname, libname) }
.Last.lib<-function(libpath) { library.dynam.unload("rbamtools",libpath) }

# #####################################################################################################################
# Generics definitions for bamReader and bamWriter
setGeneric("filename", function(object) standardGeneric("filename"))
setGeneric("isOpen",function(con,rw="") standardGeneric("isOpen"))
setGeneric("bamClose", function(object) standardGeneric("bamClose"))
setGeneric("getNextAlign",function(object) standardGeneric("getNextAlign")) # generic for bamReader and bamRange
setGeneric("bamSave",function(object,...) standardGeneric("bamSave"))       # generic for bamWriter and bamRange
setGeneric("as.data.frame",function(x,row.names=NULL,optional=FALSE,...) standardGeneric("as.data.frame")) # generic for bamRange and refSeqDict
setGeneric("as.list",function(x,...) standardGeneric("as.list"))      # Generic for conversion into list
setGeneric("getHeaderText",function(object,delim="\t")standardGeneric("getHeaderText"))  # Generic for retrieving RefData string from Objects
setGeneric("getVal",function(object,member)standardGeneric("getVal")) # Generic for Reading member from object list
setGeneric("setVal",function(object,members,values)standardGeneric("setVal")) # Generic for Writing member to object list

bamReader<-function(filename) {
  if(!file.exists(filename))
    stop("[bamReader] File '",basename(filename),"' does not exist!\n  [bamReader] Dirname: ",dirname(filename))
  return(new("bamReader",filename))
}

# #################################################################################################
# #################################################################################################
# bamReader
#
#
setClass("bamReader", representation(filename="character",reader="externalptr",index="externalptr",indexName="character"),
			validity=function(object){
				return(TRUE)
			}
)
setMethod(f="initialize", signature="bamReader",
		definition=function(.Object,filename){
			.Object@filename<-filename
			.Object@reader=.Call("bam_reader_open",path.expand(filename))
			return(.Object)
		}
)
# generic in header
setMethod("filename", "bamReader", function(object) return(object@filename))
#setGeneric("isOpen",function(con,rw="") standardGeneric("isOpen"))
setMethod("isOpen",signature="bamReader",definition=function(con,rw="") { return(!(.Call("is_nil_externalptr",con@reader,PACKAGE="rbamtools"))) })
setGeneric("index.initialized",function(object) standardGeneric("index.initialized"))
setMethod("index.initialized", signature="bamReader",definition=function(object) {return(!(.Call("is_nil_externalptr",object@index,PACKAGE="rbamtools")))})

# generic in header
setMethod(f="bamClose",signature="bamReader",definition=function(object) {
  if(!.Call("is_nil_externalptr",object@index,PACKAGE="rbamtools"))
    {.Call("bam_reader_unload_index",object@index,PACKAGE="rbamtools")}
  invisible(.Call("bam_reader_close",object@reader,PACKAGE="rbamtools"))
  })
#setGeneric("getHeaderText", function(object) standardGeneric("getHeaderText"))
setMethod(f="getHeaderText",signature="bamReader",definition=function(object) {
  .Call("bam_reader_get_header_text",object@reader,PACKAGE="rbamtools") })
setGeneric("getHeader",function(object) standardGeneric("getHeader"))
setMethod(f="getHeader",signature="bamReader",definition=function(object){
    return(new("bamHeader",.Call("bam_reader_get_header_text",object@reader,PACKAGE="rbamtools"))) })
setGeneric("getRefCount",function(object) standardGeneric("getRefCount"))
setMethod(f="getRefCount",signature="bamReader",definition=function(object) {
  return(.Call("bam_reader_get_ref_count",object@reader,PACKAGE="rbamtools"))})
setGeneric("getRefData",function(object) standardGeneric("getRefData"))
setMethod(f="getRefData",signature="bamReader",definition=function(object) {
  return(.Call("bam_reader_get_ref_data",object@reader,PACKAGE="rbamtools"))})
setGeneric("createIndex",function(object,idx_filename) standardGeneric("createIndex"))
setMethod(f="createIndex",signature="bamReader",definition=function(object,idx_filename) {
  invisible(.Call("bam_reader_create_index",object@filename,idx_filename,PACKAGE="rbamtools"))})
# setGeneric("getNextAlign",function(object) standardGeneric("getNextAlign"))
setMethod(f="getNextAlign",signature="bamReader",definition=function(object) {
  return(new("bamAlign",.Call("bam_reader_get_next_align",object@reader,PACKAGE="rbamtools"))) })
setGeneric("loadIndex",function(object,filename) standardGeneric("loadIndex"))
setMethod("loadIndex",signature="bamReader",definition=function(object,filename){
  if(!is.character(filename))
    stop("[bamReader.loadIndex] Filename must be character!\n")
  if(!file.exists(filename))
    stop("[bamReader.loadIndex] Index file \"",filename,"\" does not exist!\n")

  # Set index Slot in given bamReader object:
  # Read object name, create expression string and evaluate in parent frame
  reader<-deparse(substitute(object))
  extxt<-paste(reader,"@index<-.Call(\"bam_reader_load_index\",\"",path.expand(filename),"\",PACKAGE=\"rbamtools\")",sep="")
  eval.parent(parse(text=extxt))
  
  # Set indexName Slot in given bamReader object:
  extxt<-paste(reader,"@indexName<-'",basename(filename),"'",sep="")
  eval.parent(parse(text=extxt))
  
  # Return true if bamReader@index!=NULL (parent frame)
  extxt<-paste(".Call(\"is_nil_externalptr\",",reader,"@index,PACKAGE=\"rbamtools\")",sep="")
  invisible(!eval.parent(parse(text=extxt)))
})

setGeneric("indexName",function(object) standardGeneric("indexName"))
setMethod(f="indexName",signature="bamReader",definition=function(object) { return(object@indexName) })

setGeneric("bamSort",function(object,prefix,byName=FALSE,maxmem=100000000) standardGeneric("bamSort"))
setMethod(f="bamSort",signature="bamReader",definition=function(object,prefix="sorted",byName=FALSE,maxmem=100000000)
  {
      maxmem<-floor(maxmem)
      cat("[bamSort] Filename: ",object@filename,"\n")
      cat("[bamSort] Prefix  : ",prefix,"\n")
      cat("[bamSort] Maxmem  : ",maxmem,"\n")
      cat("[bamSort] By Name : ",byName,"\n")
      .Call("bam_reader_sort_file",object@filename,prefix,maxmem,byName,PACKAGE="rbamtools")
      cat("[bamSort] Sorting finished.\n")
      invisible(paste(prefix,".bam",sep=""))
  })


# #################################################################################################
# #################################################################################################
# bamHeader
# Description: See SAM File Format Specification (v1.4-r985) September 7,2011, Section 1.3
#
#
# #################################################################################################
# headerLine: Represents two entries: Format version (VN) and sorting order (SO)
# Valid entries for SO: unknown (default), unsorted, queryname, coordinate.

setClass("headerLine",representation(l="list"),validity=function(object){return(TRUE)})

setMethod(f="initialize",signature="headerLine",definition=function(.Object,hl="",delim="\t"){
  # Parses header line from header section
  if(!is.character(hl))
    stop("[headerLine.initialize] Argument must be string.\n")
  # hl="" or character(0)
  if((length(hl)==1 && nchar(hl)==0) || length(hl)==0)
    return(.Object)
  # Split input string into tags
  tags<-unlist(strsplit(hl,delim))
  if(tags[1]!="@HD")
    stop("[headerLine.initialize] First tag of string must be @HD!\n")
  tags<-tags[-1]
  n<-length(tags)
  tagLabs<-c("VN","SO")
  for(i in 1:n)
  {
    lab<-substr(tags[i],1,2)
    mtc<-match(lab,tagLabs)
    if(is.na(mtc))
      stop("[headerLine.initialize] Tag identifier '",lab,"' not in List!\n")
    .Object@l[[lab]]<-substr(tags[i],4,nchar(tags[i]))
  }
  return(.Object)
})

setMethod("getHeaderText",signature="headerLine",definition=function(object,delim="\t"){
  n<-length(object@l)
  if(n==0)
    return(character(0))
  
  rfstr<-character(n)
  tagnames<-names(object@l)
  for(i in 1:n)
    rfstr[i]<-paste(tagnames[i],object@l[[i]],sep=":")
  return(paste("@HD",paste(rfstr,collapse=delim),sep=delim))
})

setMethod("getVal",signature="headerLine",definition=function(object,member){
  if(!is.character(member))
    stop("[getVal.headerLine] Member must be character!\n")
  tagLabs<-c("VN","SO")
  mtc<-match(member[1],tagLabs)
  if(is.na(mtc))
    stop("[getVal.headerLine] Invalid member name!\n")
  return(object@l[[member]])
})

setMethod("setVal",signature="headerLine",definition=function(object,members,values){
  if(!is.character(members) || !is.character(values))
    stop("[setVal.headerLine] Members and values must be character!\n")
  if(length(members)!=length(values))
    stop("[setVal.headerLine] Members and values must have same length!\n")
  tagLabs<-c("VN","SO")
  mtc<-match(members,tagLabs)
  if(any(is.na(mtc)))
    stop("[setVal.headerLine] Memer names must be valid Header line entries!\n")
  n<-length(members)
  obj<-deparse(substitute(object))
  for(i in 1:n)
  {
    txt<-paste(obj,"@l$",members[i],"<-'",values[i],"'",sep="")
    eval.parent(parse(text=txt))
  }
  return(invisible())
})

setMethod("as.list",signature="headerLine",definition=function(x,...) {return(x@l)})



# #################################################################################################
# refSeqDict
# Reference Sequence Dictionary: Represents a variable number of Ref Seqs
# Valid Members (Entries for each sequence, stored in a data.frame): 
# SN Reference sequence name
# LN Reference sequence length
# AS Genome assembly identifier
# M5 MD5 checksum of the sequence
# SP Species
# UR URI of the sequence

setClass("refSeqDict",representation(df="data.frame"),
         validity=function(object) { return(TRUE)})


setMethod(f="initialize",signature="refSeqDict",definition=function(.Object,hsq="",delim="\t") {
  # Parses Reference sequence dictionary of header-text
  # Expects Vector of (internally [tab] delimited) characters, each representing one Ref-Sequence
  if(!is.character(hsq))
    stop("[refSeqDict.initialize] hsq must be character!\n")
  n<-length(hsq)
  # hsq="" or character(0)
  if((n==1 && nchar(hsq)==0) || n==0)
  {
    # Passing hsq="" creates an empty data.frame
    .Object@df<-data.frame(SN=character(0),LN=double(0),AS=character(0),M5=numeric(0),SP=character(0),UR=character(0),stringsAsFactors=F)    
    return(.Object)
  }
  .Object@df<-data.frame(SN=character(n),LN=double(n),AS=character(n),M5=numeric(n),SP=character(n),UR=character(n),stringsAsFactors=F)
  labels<-c("SN","LN","AS","M5","SP","UR")
  for(i in 1:n)
  {
    # Containes separated tags for one sequence
    seq<-((unlist(strsplit(hsq[i],delim)))[-1])
    # Contains column number in dict@df for each tag
    cols<-match(substr(seq,1,2),labels)
    m<-length(cols)
    for(j in 1:m)
    {
      txt<-substr(seq[j],4,nchar(seq[j]))
      # Empty entries are skipped (to avoid errors)
      if(nchar(txt)>0)
      {
        if(is.element(j,c(2,4)))  # 2,4 are the numeric columns of df
        {
          # Try to convert into numeric value
          numb<-suppressWarnings(as.numeric(txt))
          if(is.na(numb))
            {warning("[refSeqDict.initialize] No numeric value for \"",labels[j],"\" in \"",rsd[i],"\"!\n",sep="")}
          .Object@df[i,j]<-numb      
        }
        else
        { .Object@df[i,j]<-txt }       
      }
    }
  }
  return(.Object)
})

setMethod(f= "[",signature="refSeqDict",definition=function(x,i,j,drop){
  return(x@df[i,j])
})
setReplaceMethod(f="[",signature="refSeqDict",definition=function(x,i,j,value){
  x@df[i,j]<-value
  if(!is.numeric(x@df[,2]))
  {
    warning("[refSeqDict ReplaceMethod\"[\"] Non-numeric entry in Column 2 (LN=Sequence Length)!\n")
    x@df[,2]<-trunc(as.numeric(x@df[,2]))
  }
  if(!is.numeric(x@df[,4]))
  {
    warning("[refSeqDict ReplaceMethod\"[\"] Non-numeric entry in Column 4 (M5=MD5 checksum)!\n")
    x@df[,4]<-trunc(as.numeric(x@df[,4]))
  }
  return(x)
})

setMethod(f="dim",signature="refSeqDict",definition=function(x){return(dim(x@df))})
setMethod("as.data.frame",signature="refSeqDict",definition=function(x,row.names=NULL,optional=FALSE,...)
  { return(as.data.frame(x@df,row.names=row.names(x@df),optional=optional,...)) })

setGeneric("removeRows",function(object,rows)standardGeneric("removeRows"))
setMethod("removeRows",signature="refSeqDict",definition=function(object,rows){
  # Removes given rows (=Sequences) from Dictionary so they are excluded from header
  dfRows<-1:(dim(object@df)[1])
  rmRows<-sort(unique(rows))
  if(any(is.na(match(rmRows,dfRows))))
    stop("[removeRows] Trying to remove rows not present in data!\n")
  
  # Construct command string for eval in parent.frame
  rmv<-paste("c(",paste(rmRows,collapse=","),")",sep="")
  dictdf<-paste(deparse(substitute(object)),"@df",sep="")
  eval.parent(parse(text=paste(dictdf,"<-",dictdf,"[-",rmv,",]",sep="")))
  return(invisible())
})

setGeneric("addSeq",function(object,SN,LN,AS="",M5=0,SP="",UR="")standardGeneric("addSeq"))
setMethod("addSeq",signature="refSeqDict",definition=function(object,SN,LN,AS="",M5=0,SP="",UR=""){
  # Appends new Sequence (row) at the end
  newRow<-paste("data.frame(SN='",SN,"',LN=",LN,",AS='",AS,"',M5=",M5,",SP='",SP,"',UR='",UR,"')",sep="")
  dictdf<-paste(deparse(substitute(object)),"@df",sep="")
  eval.parent(parse(text=paste(dictdf,"<-rbind(",dictdf,",",newRow,")",sep="")))
  return(invisible())
})

#setGeneric("getHeaderText",function(object)standardGeneric("getHeaderText"))
setMethod("getHeaderText",signature="refSeqDict",definition=function(object,delim="\t"){ 
  # Returns Ref Data String (can be used for creating new BAM file via bamWriter)
  labels<-c("SN","LN","AS","M5","SP","UR")
  n<-dim(object@df)[1]
  if(n==0)
    return(character(0))
  
  seqs<-character(n)
  nCols<-6
  seqEntries<-character(nCols)
  for(i in 1:n)
  {
    # Find out which columns must be present in output for this seq
    colsInUse<-rep(FALSE,nCols)
    for(j in c(2,4))
    {
      if(object@df[i,j]>0)
        colsInUse[j]<-TRUE      
    }
    for(j in c(1,3,5,6))
    {
      if(nchar(object@df[i,j])>0)
        colsInUse[j]<-TRUE
    }
    # First two columns (SN,LN) must always be present!
    if(!colsInUse[1]||!colsInUse[2])
      stop("[getHeaderText.refSeqDict] Mandantory entry in column(s)",
           paste(labels[!(colsInUse[1:2])],collapse=",")," missing!\n")
    
    #Write output for each used column
    for(j in 1:nCols)
      seqEntries[j]<-paste(labels[j],":",object@df[i,j],sep="")
    seqs[i]<-paste("@SQ",paste(seqEntries[colsInUse],collapse=delim),sep=delim)
  }
  return(paste(seqs,collapse="\n"))
})

# Return first or last part of refSeqDict data.frame
# S3 Generic is supplied via importFrom in NAMESPACE
setGeneric("head",function(x,...) standardGeneric("head"))
setMethod("head","refSeqDict",function(x,n=6L,...) {
  stopifnot(length(n) == 1L)
  if (n < 0L)
    n<-0
  else
    n<-min(n, dim(x@df)[1])
  if(n == 0L)
    return(x@df[integer(0),])
  else
    return(x@df[1:n,])
})
# S3 Generic is supplied via imporFrom in NAMESPACE
setGeneric("tail",function(x,...) standardGeneric("tail"))
setMethod("tail","refSeqDict",definition=function(x,n=6L,...) {
  stopifnot(length(n) == 1L)
  nRows<-dim(x@df)[1]
  if (n < 0L)
    n<-0
  else
    n<-min(n,nRows)
  if(n == 0L)
    return(x@df[integer(0),])
  else
  {
    n<-nRows-n+1
    return(x@df[nRows:n,])
  }
})


###################################################################################################
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
# PL Sequencing Platform (CAPILLARY,LS454,ILLUMINA,SOLID,HELICOS,IONTORRENT or PACBIO)
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

# setGeneric("getHeaderText",function(object)standardGeneric("getHeaderText"))
setMethod("getHeaderText",signature="headerReadGroup",definition=function(object,delim="\t") {
  n<-length(object@l)
  if(n==0)
    return(character(0))
  rfstr<-character(n)
  for(i in 1:n)
    rfstr[i]<-paste(names(object@l)[i],object@l[[i]],sep=":")
  return(paste("@RG",paste(rfstr,collapse=delim),sep=delim))
})

# setGeneric("getVal",function(object,member)standardGeneric("getVal"))
setMethod("getVal",signature="headerReadGroup",definition=function(object,member) {
  if(!is.character(member))
    stop("[getVal.headerReadGroup] Member must be character!\n")
  tagLabs<-c("ID","CN","DS","DT","FO","KS","LB","PG","PI","PL","PU","SM")
  mtc<-match(member[1],tagLabs)
  if(is.na(mtc))
    stop("[getVal.headerReadGroup] Invalid member name!\n")
  return(object@l[[member]])
})

#setGeneric("setVal",function(object,members,values)standardGeneric("setVal"))
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

###################################################################################################
# headerProgram



setClass("headerProgram",representation(l="list"),validity=function(object){return(TRUE)})

setMethod(f="initialize",signature="headerProgram", definition=function(.Object,hp="",delim="\t"){
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

# setGeneric("getHeaderText",function(object)standardGeneric("getHeaderText"))
setMethod("getHeaderText",signature="headerProgram",definition=function(object,delim="\t") {
  n<-length(object@l)
  if(n==0)
    return(character(0))
  
  rfstr<-character(n)
  for(i in 1:n)
    rfstr[i]<-paste(names(object@l)[i],object@l[[i]],sep=":")
  return(paste("@PG",paste(rfstr,collapse=delim),sep=delim))
})

# setGeneric("getVal",function(object,member)standardGeneric("getVal"))
setMethod("getVal",signature="headerProgram",definition=function(object,member) {
  if(!is.character(member))
    stop("[getVal.headerProgram] Member must be character!\n")
  tagLabs<-c("ID","PN","CL","PP","VN")
  mtc<-match(member[1],tagLabs)
  if(is.na(mtc))
    stop("[getVal.headerProgram] Invalid member name!\n")
  return(object@l[[member]])
})

#setGeneric("setVal",function(object,members,values)standardGeneric("setVal"))
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

# #################################################################################################
# bamHeader: Container for complete BAM Header segment

setClass("bamHeader",representation(head="headerLine",dict="refSeqDict",
                    group="headerReadGroup",prog="headerProgram",com="character"))
setMethod(f="initialize",signature="bamHeader", definition=function(.Object,bh="",delim="\n"){
  # Parses Header data (as reported by getHeaderText)
  # See Sam Format Specificatioin 1.3 (Header Section)
  if(!is.character(bh))
    stop("[bamHeader.initialize] Argument must be string.\n")
  
  # Create empty header Set (so it's legal to call getHeaderText')
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

setGeneric("getHeaderLine",function(object) standardGeneric("getHeaderLine"))
setMethod(f="getHeaderLine",signature="bamHeader",definition=function(object) {return(object@head)})
setGeneric("getRefSeqDict",function(object) standardGeneric("getRefSeqDict"))
setMethod(f="getRefSeqDict",signature="bamHeader",definition=function(object) {return(object@dict)})
setGeneric("getHeaderReadGroup",function(object)standardGeneric("getHeaderReadGroup"))
setMethod(f="getHeaderReadGroup",signature="bamHeader",definition=function(object){return(object@group)})
setGeneric("getHeaderProgram",function(object)standardGeneric("getHeaderProgram"))
setMethod(f="getHeaderProgram",signature="bamHeader",definition=function(object){return(object@prog)})

setMethod("getHeaderText",signature="bamHeader",definition=function(object,delim="\n") {
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
    cm<-paste(paste("HD",object@com,sep="\t"),collapse=delim)
  else
    cm<-character(0)
  return(paste(hd,dt,gp,pg,cm,sep=""))
})


# #################################################################################################
# #################################################################################################
# bamWriter class
# Encapsulates an write-opened Connection to a BAM-file.
#

setClass("bamWriter", representation(filename="character",writer="externalptr"),
			validity=function(object){
				return(TRUE)
			}
)

setMethod(f="initialize", signature="bamWriter",
		definition=function(.Object,reader,filename){
			.Object@filename<-filename
      .Object@writer<-.Call("bam_writer_open",reader@reader,filename,PACKAGE="rbamtools")
			return(.Object)
		}
)

bamWriter<-function(reader,filename) { return(new("bamWriter",reader,filename)) }

# generic in header
setMethod("filename", "bamWriter", function(object) return(object@filename))
#setGeneric("isOpen",function(con,rw="") standardGeneric("isOpen"))
setMethod("isOpen",signature="bamWriter",definition=function(con,rw="") { return(!(.Call("is_nil_externalptr",con@writer,PACKAGE="rbamtools"))) })


# generic in header
setMethod(f="bamClose",signature="bamWriter",definition=function(object) { invisible(.Call("bam_writer_close",object@writer,PACKAGE="rbamtools"))})
# setGeneric("bamSave",function(object,...) standardGeneric("bamSave"))
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

# #################################################################################################
# #################################################################################################
# bamRange
#
# Encapsulates a bunch of Alignment datasets that typically have been read from a defined region
# on the reference genome. Technically, the alignments are stored in a (C-implemented) double
# linked list.
# bamRange objects can be created by a reading procedure on an indexed BAM-file. The alignments can
# iterated, readed, written, deleted and added. bamRange objects can be written to a BAM-file via
# an Instance of bamWriter.

bamRange<-function(reader,coords) { return(new("bamRange",reader,coords))}

setClass("bamRange",representation(range="externalptr"),
         validity=function(object) { return(ifelse(is.null(.Object@range),FALSE,TRUE)) })

setMethod(f="initialize",signature="bamRange",
          definition=function(.Object,reader,coords) {
              if(!is(reader,"bamReader"))
                stop("bamRange initialize: reader must be an instance of bamReader!\n")
              if(length(coords)!=3)
                stop("bamRange initialize: coords must be 3-dim numeric (ref,start,stop)!\n")  
              if(is.null(reader@index))
                stop("bamRange initialize: bamReader must have initialized index!\n")
              .Object@range<-.Call("bam_range_fetch",reader@reader,reader@index,trunc(coords),PACKAGE="rbamtools")
              return(.Object)
              })
#setGeneric("as.data.frame",function(x,row.names=NULL,optional=FALSE,...) standardGeneric("as.data.frame"))
setMethod("as.data.frame",signature="bamRange",definition=function(x,row.names=NULL,optional=FALSE,...) {
        df<-.Call("bam_range_get_align_df",x@range,PACKAGE="rbamtools")
        as.data.frame(df,row.names=row.names(df),optional=optional,...)
        })

setGeneric("size",function(object) standardGeneric("size"))
setMethod("size",signature="bamRange",definition=function(object) { .Call("bam_range_get_size",object@range,PACKAGE="rbamtools")} )

#setGeneric("getNextAlign",function(object) standardGeneric("getNextAlign"))
setMethod("getNextAlign",signature="bamRange",definition=function(object) { return(new("bamAlign",.Call("bam_range_get_next_align",object@range,PACKAGE="rbamtools")))})

setGeneric("getPrevAlign",function(object) standardGeneric("getPrevAlign"))
setMethod("getPrevAlign",signature="bamRange",definition=function(object) { return(new("bamAlign",.Call("bam_range_get_prev_align",object@range,PACKAGE="rbamtools")))})

setGeneric("windBack", function(object) standardGeneric("windBack"))
setMethod("windBack",signature="bamRange",definition=function(object) {invisible(.Call("bam_range_wind_back",object@range,PACKAGE="rbamtools"))})

setGeneric("push_back",function(object,value) standardGeneric("push_back"))
setMethod("push_back",signature="bamRange",definition=function(object,value)
  {
    if(!is(value,"bamAlign"))
      stop("push_back.bamRange: pushed object must be of class \"bamAlign\"\n")
    .Call("bam_range_push_back",object@range,value@align,PACKAGE="rbamtools")
  })

setGeneric("pop_back",function(object) standardGeneric("pop_back"))
setMethod("pop_back",signature="bamRange",definition=function(object){.Call("bam_range_pop_back",object@range,PACKAGE="rbamtools") })

setGeneric("push_front",function(object,value) standardGeneric("push_front"))
setMethod("push_front",signature="bamRange",definition=function(object,value)
  {
    if(!is(value,"bamAlign"))
      stop("push_front.bamRange: pushed object must be of class \"bamAlign\"\n")
    .Call("bam_range_push_front",object@range,value@align,PACKAGE="rbamtools")
  })

setGeneric("pop_front",function(object) standardGeneric("pop_front"))
setMethod("pop_front",signature="bamRange",definition=function(object) {.Call("bam_range_pop_front",object@range,PACKAGE="rbamtools")})

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


# #################################################################################################
# #################################################################################################
# bamAlign
#
# bamAlign encapsulates all contained data in a single dataset in a BAM-file. bamAlign objects can
# be read from a bamReader instance and written to a bamWriter instance. All contained data can be
# read and written via accessor functions.
#

setClass("bamAlign", representation(align="externalptr"),
			validity=function(object){return(ifelse(is.null(.Object@align,FALSE,TRUE)))})

setMethod(f="initialize", signature="bamAlign",
		definition=function(.Object,align=NULL){
			.Object@align<-align
			return(.Object)
		}
)

# bamAlign Member Reader functions
setGeneric("name",function(object) standardGeneric("name"))
setMethod(f="name",signature="bamAlign",definition=function(object) { .Call("bam_align_get_name",object@align,PACKAGE="rbamtools") })
setGeneric("refID",function(object) standardGeneric("refID"))
setMethod(f="refID",signature="bamAlign",definition=function(object) {return(.Call("bam_align_get_refid",object@align,PACKAGE="rbamtools"))})
setGeneric("position",function(object) standardGeneric("position"))
setMethod(f="position",signature="bamAlign",definition=function(object) {return(.Call("bam_align_get_position",object@align,PACKAGE="rbamtools"))})
setGeneric("cigarData",function(object) standardGeneric("cigarData"))
setMethod(f="cigarData",signature="bamAlign",definition=function(object) { .Call("bam_align_get_cigar_df",object@align,PACKAGE="rbamtools")})
setGeneric("mateRefID",function(object) standardGeneric("mateRefID"))
setMethod(f="mateRefID",signature="bamAlign",definition=function(object) { .Call("bam_align_get_mate_refid",object@align,PACKAGE="rbamtools")})
setGeneric("matePosition",function(object) standardGeneric("matePosition"))
setMethod(f="matePosition",signature="bamAlign",definition=function(object) { .Call("bam_align_get_mate_position",object@align,PACKAGE="rbamtools")})
setGeneric("insertSize",function(object) standardGeneric("insertSize"))
setMethod(f="insertSize",signature="bamAlign",definition=function(object) { .Call("bam_align_get_insert_size",object@align,PACKAGE="rbamtools")})
setGeneric("mapQuality",function(object) standardGeneric("mapQuality"))
setMethod(f="mapQuality",signature="bamAlign",definition=function(object) { .Call("bam_align_get_map_quality",object@align,PACKAGE="rbamtools")})
setGeneric("readBases",function(object) standardGeneric("readBases"))
setMethod(f="readBases",signature="bamAlign",definition=function(object) {return(.Call("bam_align_get_read_bases",object@align,PACKAGE="rbamtools"))})
setGeneric("qualities",function(object) standardGeneric("qualities"))
setMethod(f="qualities",signature="bamAlign",definition=function(object) {return(.Call("bam_align_get_qualities",object@align,PACKAGE="rbamtools"))})

# ##############################################################
# Queries against alignment flag (Readers and Accessors)

setGeneric("pcrORopt_duplicate", function(object) standardGeneric("pcrORopt_duplicate"))
setGeneric("pcrORopt_duplicate<-", function(object,value) standardGeneric("pcrORopt_duplicate<-"))
setMethod("pcrORopt_duplicate", "bamAlign", function(object) return(.Call("bam_align_is_pcr_or_optical_dup",object@align,PACKAGE="rbamtools")))
setReplaceMethod(f="pcrORopt_duplicate", signature="bamAlign",
   definition=function(object,value){
       if(!is.logical(value))
           stop("class bamReader, Duplicate setter: value must be boolean")
         .Call("bam_align_set_is_pcr_or_optical_dup",object@align,value,PACKAGE="rbamtools")
       return(object)
   }
)

setGeneric("failedQC", function(object) standardGeneric("failedQC"))
setGeneric("failedQC<-", function(object,value) standardGeneric("failedQC<-"))
setMethod("failedQC", "bamAlign", function(object) return(.Call("bam_align_fail_qc",object@align,PACKAGE="rbamtools")))
setReplaceMethod(f="failedQC", signature="bamAlign",
   definition=function(object,value){
       if(!is.logical(value))
           stop("class bamReader, failedQC setter: value must be boolean")
       .Call("bam_align_set_fail_qc",object@align,value,PACKAGE="rbamtools")
       return(object)
   }
)

setGeneric("firstInPair", function(object) standardGeneric("firstInPair"))
setGeneric("firstInPair<-", function(object,value) standardGeneric("firstInPair<-"))
setMethod("firstInPair", "bamAlign", function(object) return(.Call("bam_align_is_first_in_pair",object@align,PACKAGE="rbamtools")))
setReplaceMethod(f="firstInPair", signature="bamAlign",
   definition=function(object,value){
       if(!is.logical(value))
           stop("class bamReader, FirstInPair setter: value must be boolean")
       .Call("bam_align_set_is_first_in_pair",object@align,value,PACKAGE="rbamtools")
       return(object)
   }
)

setGeneric("secondInPair", function(object) standardGeneric("secondInPair"))
setGeneric("secondInPair<-", function(object,value) standardGeneric("secondInPair<-"))
setMethod("secondInPair", "bamAlign", function(object) return(.Call("bam_align_is_second_in_pair",object@align,PACKAGE="rbamtools")))
setReplaceMethod(f="secondInPair", signature="bamAlign",
   definition=function(object,value){
       if(!is.logical(value))
           stop("class bamReader, secondInPair setter: value must be boolean")
       .Call("bam_align_set_is_second_in_pair",object@align,value,PACKAGE="rbamtools")
       return(object)
   }
)

setGeneric("unmapped", function(object) standardGeneric("unmapped"))
setGeneric("unmapped<-", function(object,value) standardGeneric("unmapped<-"))
setMethod("unmapped", "bamAlign", function(object) return(.Call("bam_align_is_unmapped",object@align,PACKAGE="rbamtools")))
setReplaceMethod(f="unmapped", signature="bamAlign",
   definition=function(object,value){
       if(!is.logical(value))
           stop("class bamReader, unmapped setter: value must be boolean")
       .Call("bam_align_set_is_unmapped",object@align,value,PACKAGE="rbamtools")
       return(object)
   }
)

setGeneric("mateUnmapped", function(object) standardGeneric("mateUnmapped"))
setGeneric("mateUnmapped<-", function(object,value) standardGeneric("mateUnmapped<-"))
setMethod("mateUnmapped", "bamAlign", function(object) return(.Call("bam_align_mate_is_unmapped",object@align,PACKAGE="rbamtools")))
setReplaceMethod(f="mateUnmapped", signature="bamAlign",
   definition=function(object,value){
       if(!is.logical(value))
           stop("class bamReader, mateUnmapped setter: value must be boolean")
       .Call("bam_align_set_mate_is_unmapped",object@align,value,PACKAGE="rbamtools")
       return(object)
   }
)

setGeneric("reverseStrand", function(object) standardGeneric("reverseStrand"))
setGeneric("reverseStrand<-", function(object,value) standardGeneric("reverseStrand<-"))
setMethod("reverseStrand", "bamAlign", function(object) return(.Call("bam_align_strand_reverse",object@align,PACKAGE="rbamtools")))
setReplaceMethod(f="reverseStrand", signature="bamAlign",
   definition=function(object,value){
       if(!is.logical(value))
           stop("class bamReader, reverseStrand setter: value must be boolean")
       .Call("bam_align_set_strand_reverse",object@align,value,PACKAGE="rbamtools")
       return(object)
   }
)

setGeneric("mateReverseStrand", function(object) standardGeneric("mateReverseStrand"))
setGeneric("mateReverseStrand<-", function(object,value) standardGeneric("mateReverseStrand<-"))
setMethod("mateReverseStrand", "bamAlign", function(object) return(.Call("bam_align_mate_strand_reverse",object@align,PACKAGE="rbamtools")))
setReplaceMethod(f="mateReverseStrand", signature="bamAlign",
   definition=function(object,value){
       if(!is.logical(value))
           stop("class bamReader, mateReverseStrand setter: value must be boolean")
       .Call("bam_align_set_mate_strand_reverse",object@align,value,PACKAGE="rbamtools")
       return(object)
   }
)

setGeneric("paired", function(object) standardGeneric("paired"))
setGeneric("paired<-", function(object,value) standardGeneric("paired<-"))
setMethod("paired", "bamAlign", function(object) return(.Call("bam_align_is_paired",object@align,PACKAGE="rbamtools")))
setReplaceMethod(f="paired", signature="bamAlign",
   definition=function(object,value){
       if(!is.logical(value))
           stop("class bamReader, paired setter: value must be boolean")
       .Call("bam_align_set_is_paired",object@align,value,PACKAGE="rbamtools")
       return(object)
   }
)


setGeneric("properPair", function(object) standardGeneric("properPair"))
setGeneric("properPair<-", function(object,value) standardGeneric("properPair<-"))
setMethod("properPair", "bamAlign", function(object) return(.Call("bam_align_mapped_in_proper_pair",object@align,PACKAGE="rbamtools")))
setReplaceMethod(f="properPair", signature="bamAlign",
   definition=function(object,value){
       if(!is.logical(value))
           stop("class bamReader, properPair setter: value must be boolean")
       .Call("bam_align_set_mapped_in_proper_pair",object@align,value,PACKAGE="rbamtools")
       return(object)
   }
)

setGeneric("secondaryAlign", function(object) standardGeneric("secondaryAlign"))
setGeneric("secondaryAlign<-", function(object,value) standardGeneric("secondaryAlign<-"))
setMethod("secondaryAlign", "bamAlign", function(object) return(.Call("bam_align_is_secondary_align",object@align,PACKAGE="rbamtools")))
setReplaceMethod(f="secondaryAlign", signature="bamAlign",
   definition=function(object,value){
       if(!is.logical(value))
           stop("class bamReader, SecondaryAlign setter: value must be boolean")
       .Call("bam_align_set_is_secondary_align",object@align,value,PACKAGE="rbamtools")
       return(object)
   }
)

setGeneric("flag", function(object) standardGeneric("flag"))
setGeneric("flag<-", function(object,value) standardGeneric("flag<-"))
setMethod("flag", "bamAlign", function(object) return(.Call("bam_align_get_flag",object@align,PACKAGE="rbamtools")))
setReplaceMethod(f="flag", signature="bamAlign",
   definition=function(object,value){
       if(!is.integer(value))
           stop("class bamReader, flag setter: value must be boolean")
       .Call("bam_align_set_flag",object@align,value,PACKAGE="rbamtools")
       return(object)
   }
)
