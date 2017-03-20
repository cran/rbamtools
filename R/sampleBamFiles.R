

setGeneric("extractGeneRegions",
           function(src, trg, gl) standardGeneric("extractGeneRegions"))

# + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + #
# CLASS sampleBamFiles
# + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + #

.sampleBamFiles <- setClass("sampleBamFiles",
    representation(bamFiles="character",
        bamIdxFiles="character",
        nAligns="numeric",
        group="factor",
        label="character",
        length="integer",
        ev="environment"),
    prototype=prototype(
        bamFiles="",
        bamIdxFiles="",
        nAligns=1,
        group=factor(),
        length=0L,
        ev=new.env()),
    validity=function(object){
        if(length(object@bamFiles) != length(object@group))
            return("bamFiles and group must have equal length.")

        if(length(object@bamFiles) != length(object@bamIdxFiles))
            return("bamFiles and bamIdxFiles must have equal length")
    }
)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# bamFiles      : BAM file locations
# bamIdxFiles   : BAM index file locations
# label         : (short) Sample label. Used for printing and plotting.
# group         : Group assignments for BAM files
# nAligns       : Total aligned reads for each sample
#                 (used for align-depth normalization over multiple samples)
# length        : Vector length for bamFiles, bamIdxFiles, label, group,
#                   nAligns. Shall not be changed after object creation
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #

setMethod("initialize", "sampleBamFiles", function(.Object)
{
    .Object@ev=new.env()
    return(.Object)
})


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# length
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
setMethod("length", "sampleBamFiles", function(x){
    return(length(x@bamFiles))
})

setMethod("[", signature="sampleBamFiles", function(x, i)
{
    i <- sort(unique(as.integer(i)))

    if(i[length(i)] > length(x))
        stop("Out of bounds index (> length(x))!")

    res <- sampleBamFiles(length(i))
    res@bamFiles <- x@bamFiles[i]
    res@bamIdxFiles <- x@bamIdxFiles[i]
    res@nAligns <- x@nAligns[i]
    res@group <- factor(x@group[i])
    res@label <- x@label[i]
    res@length <- length(i)

    if(exists("groupTable", envir=x@ev))
        assign("groupTable", x@ev$groupTable[i, ], envir=res@ev)

    return(res)
})

setMethod("show", "sampleBamFiles", function(object)
{
    cat("An object of class \"", class(object), "\"\n", sep="")
    cat("Number of Files :", length(object), "\n")
    cat("Groups          :", levels(object@group), "\n")
    cat("Group table:")
    print(table(object@group))
})

setGeneric("sampleBamFiles", function(object) standardGeneric("sampleBamFiles"))

setMethod("sampleBamFiles", "integer", function(object){
    len <- object[1]
    res <- .sampleBamFiles()
    res@length <- len
    res@bamFiles <- character(len)
    res@bamIdxFiles <- character(len)
    res@label <- character(len)
    res@nAligns <- numeric(len)
    return(res)
})

setMethod("sampleBamFiles", "numeric", function(object){
    return(sampleBamFiles(as.integer(object)))
})

setMethod("sampleBamFiles", "character", function(object){
    bs <- sampleBamFiles(length(object))
    bamFiles(bs) <- object
    return(bs)
})

setMethod("sampleBamFiles", "factor", function(object)
        { return(sampleBamFiles(as.character(object)))})

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# bamFiles
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
setGeneric("bamFiles", function(object) standardGeneric("bamFiles"))
setMethod("bamFiles", "sampleBamFiles", function(object){
    return(object@bamFiles)
})

setGeneric("bamFiles<-",
        function(object, value) standardGeneric("bamFiles<-"))

setReplaceMethod("bamFiles", c("sampleBamFiles", "character"),
    function(object, value)
    {
        if(length(value) != length(object))
            stop("length of value and replacement must be equal!")

        object@bamFiles <- value

        # Only fill when no idx file name given
        if(nchar(object@bamIdxFiles[1])==0)
            object@bamIdxFiles <- paste(value, "bai", sep=".")

        return(object)
})



# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# bamIdxFiles
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
setGeneric("bamIdxFiles", function(object) standardGeneric("bamIdxFiles"))
setMethod("bamIdxFiles", "sampleBamFiles", function(object){
    return(object@bamIdxFiles)
})

setGeneric("bamIdxFiles<-",
           function(object, value) standardGeneric("bamIdxFiles<-"))

setReplaceMethod("bamIdxFiles", c("sampleBamFiles", "character"),
    function(object, value)
    {
        if(length(value) != length(object))
            stop("length of value and replacement must be equal!")

        object@bamIdxFiles <- value
        return(object)
    }
)



# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# nAligns
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
setMethod("nAligns", "sampleBamFiles", function(object){
    return(object@nAligns)
})

setGeneric("nAligns<-", function(object, value) standardGeneric("nAligns<-"))
setReplaceMethod("nAligns", c("sampleBamFiles","numeric"), function(object, value){
    if(length(value) != length(object))
        stop("length of value and replacement must be equal!")

    object@nAligns <- as.numeric(value)
    return(object)
})

setReplaceMethod("nAligns", c("sampleBamFiles", "integer"), function(object, value){
    nAligns(object) <- as.numeric(value)
    return(object)
})


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# sampleLabels
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
setGeneric("sampleLabels", function(object) standardGeneric("sampleLabels"))
setMethod("sampleLabels", "sampleBamFiles", function(object){
    return(object@label)
})

setGeneric("sampleLabels<-",
                    function(object, value) standardGeneric("sampleLabels<-"))
setReplaceMethod("sampleLabels", c("sampleBamFiles", "character"),
        function(object, value)
{
    if(length(value) != length(object))
        stop("length of value and replacement must be equal!")

    object@label <- value
    return(object)
})

setReplaceMethod("sampleLabels", c("sampleBamFiles", "factor"),
        function(object, value)
{
    sampleLabels(object) <- as.character(value)
    return(object)
})

setReplaceMethod("sampleLabels", c("sampleBamFiles", "integer"),
                 function(object, value)
                 {
                     sampleLabels(object) <- as.character(value)
                     return(object)
                 })

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# sampleGroups
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
setGeneric("sampleGroups", function(object) standardGeneric("sampleGroups"))
setMethod("sampleGroups", "sampleBamFiles", function(object){
    return(object@group)
})

setGeneric("sampleGroups<-", function(object, value)
                standardGeneric("sampleGroups<-"))

setReplaceMethod("sampleGroups", c("sampleBamFiles", "factor"),
    function(object, value)
    {
        if(length(value) != length(object))
            stop("length of value and replacement must be equal!")
        object@group <- value
        return(object)
    }
)


setReplaceMethod("sampleGroups", c("sampleBamFiles", "character"),
    function(object, value)
    {
        sampleGroups(object) <- factor(value)
        return(object)
    }
)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# groupTable
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
setGeneric("groupTable", function(object) standardGeneric("groupTable"))
setMethod("groupTable", "sampleBamFiles", function(object){
    if(exists("groupTable", envir=object@ev))
        return(object@ev$groupTable)
    return(NULL)
})

setGeneric("groupTable<-",
                function(object, value) standardGeneric("groupTable<-"))

setReplaceMethod("groupTable", c("sampleBamFiles", "data.frame"),
    function(object, value)
    {
        if(nrow(value)!=length(object))
            warning("nrow(value) should be equal to length(object)!")

        object@ev$groupTable <- value
        return(object)
    }
)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# checkBamFiles
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #

setGeneric("checkBamFiles", function(x) standardGeneric("checkBamFiles"))

setMethod("checkBamFiles", "sampleBamFiles", function(x)
{
    # - - - - - - - - - - - - - - - - - - - - - - - - - - #
    # Check File existence
    # - - - - - - - - - - - - - - - - - - - - - - - - - - #
    if(!all(file.exists(x@bamFiles))){
        cat("[checkBamFiles] BAM files not found!\n")
        return(FALSE)
    }

    if(!all(file.exists(x@bamIdxFiles))){
        cat("[checkBamFiles] BAM Index files not found!\n")
        return(FALSE)
    }
    cat("[checkBamFiles] File check ok")

    # - - - - - - - - - - - - - - - - - - - - - - - - - - #
    # Check File content (SN: sequence names in RefData)
    # - - - - - - - - - - - - - - - - - - - - - - - - - - #
    frd <- bamReader(x@bamFiles[1], x@bamIdxFiles[1])
    ffd <- getRefData(frd)

    n <- length(x)
    nfok <- 1
    nrok <- 1
    if(n > 1)
    {
        for(i in 2:n)
        {
            rd <- bamReader(x@bamFiles[i], x@bamIdxFiles[i])
            if(!isOpen(rd))
            {
                cat("\n[checkBamFiles] Cannot open file No",
                        i, ".\n")
                ok <- FALSE
                section_ok <- FALSE
            }else{
                nfok <- nfok + 1

                fd <- getRefData(rd)
                if(all(!is.na(match(fd$SN, fd$SN))))
                {
                    cat("\r[checkBamFiles] BAM file",
                            format(i, w=2), " RefData OK.")
                    nrok <- nrok + 1
                }else{
                    cat("\n[checkBamFiles] Missing matches for SN (RefData) in File", i, "\n")
                    ok <- FALSE
                    section_ok <- FALSE
                }
            }
        }
    }
    cat("\n")
    cat("[checkBamFiles]", nfok, " BAM Files OK.\n")
    cat("[checkBamFiles]", nrok, " RefData   OK.\n")
    cat("[checkBamFiles] BAM + RefData content OK.\n")

    # - - - - - - - - - - - - - - - - - - - - - - - - - - #
    # Check complete
    # - - - - - - - - - - - - - - - - - - - - - - - - - - #
    cat("[checkBamFiles] Check OK.\n")
    return(invisible(TRUE))
})


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# bamCountAll
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
setMethod("bamCountAll", c("sampleBamFiles", "logical"),
    function(object, verbose=FALSE)
    {
        n <- length(object)

        if(!all(file.exists(object@bamFiles)))
            stop("BAM file not found!")

        if(!all(file.exists(object@bamIdxFiles)))
            stop("BAM index file not found!")

        reader <- bamReader(object@bamFiles[1], object@bamIdxFiles[1])
        bc <- bamCountAll(reader, verbose=verbose)
        bamClose(reader)
        bc$file <- 1

        if(n > 1)
        {
            for(i in 2:n)
            {
                reader <- bamReader(object@bamFiles[i], object@bamIdxFiles[i])
                count <- bamCountAll(reader, verbose=verbose)
                bamClose(reader)
                count$file <- i

                bc <- rbind(bc, count)
            }
        }

        object@ev$bamCounts <- bc
        return(as.numeric(tapply(bc$nAligns, bc$file, sum, na.rm=TRUE)))
    }
)


setMethod("bamCountAll", c("sampleBamFiles", "missing"), function(object, verbose){
    return(bamCountAll(object, FALSE))
})


# + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + #
# CLASS geneAlignDepth
# + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + #

# length: ncol(ald)

.geneAlignDepth <- setClass("geneAlignDepth",
    representation(ald="matrix",
        gapSites="data.frame",
        gene_id="character",
        gene_name="character",
        seq_name="character",
        strand="character",
        coords="integer"),
    contains="sampleBamFiles"
)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# length
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
setMethod("length", "geneAlignDepth", function(x){
    callNextMethod(x)
})

setMethod("[", signature="geneAlignDepth", function(x, i)
{
    i <- sort(unique(as.integer(i)))

    res <- callNextMethod(x, i)
    res@ald <- x@ald[i, ]

    res@gene_id <- x@gene_id
    res@gene_name <- x@gene_name
    res@seq_name <- x@seq_name
    res@strand <- x@strand

    if(exists("groupTable", envir=x@ev))
        assign("groupTable", x@ev$groupTable[i, ], envir=res@ev)

    return(res)
})




setMethod("initialize", "geneAlignDepth", function(.Object, basa)
{
    if(!is(basa, "sampleBamFiles"))
        stop("sampleBamFiles object needed")

    .Object@ald <- matrix()
    .Object@group <- basa@group
    .Object@bamFiles <- basa@bamFiles
    .Object@bamIdxFiles <- basa@bamIdxFiles

    if(any(basa@nAligns==0))
    {
        warning("0 nAligns present. Set nAligns to 1")
        .Object@nAligns <- rep(1, length(basa@bamFiles))
    }else{
        .Object@nAligns <- basa@nAligns
    }

    .Object@ev=new.env()
    if(exists("exp", envir=basa@ev))
        assign("exp", basa@ev$exp, envir = .Object@ev)

    return(.Object)
})



setGeneric("geneAlignDepth", function(basa, gm, check=FALSE)
                                        standardGeneric("geneAlignDepth"))

setMethod("geneAlignDepth", c("sampleBamFiles", "geneModel"),
    function(basa, gm, check=FALSE)
    {
        if(check[1])
        {
            if(!checkBamFiles(basa))
                stop("Check BAM files failed!")
        }

        res <- .geneAlignDepth(basa)

        # Copy member data from incoming geneModel
        res@gene_id <- gm@gene_id
        res@gene_name <- gm@gene_name
        res@seq_name <- gm@seq_name
        res@strand <- gm@strand
        res@coords <- gm@coords

        # Construct alignment depth matrix
        # Rows: Reference sequence position
        # Columns: Samples
        res@ald <- matrix(0L, nrow=res@coords[2] - res@coords[1] + 1,
                          ncol=length(basa))

        colnames(res@ald) <- 1:length(basa)
        rownames(res@ald) <- (res@coords[1]:res@coords[2])
        nFiles <- length(res@bamFiles)

        # Create empty gapSiteList objects for each group for later merging
        groups <- as.numeric(basa@group)
        nGroups <- length(levels(basa@group))
        gsl <- lapply(1:nGroups, siteList)

        # One global gapSiteList
        glsl <- siteList()


        # Fill alignment depth matrix
        for(i in 1:nFiles)
        {
            reader <- bamReader(res@bamFiles[i], res@bamIdxFiles[i])
            rd <- getRefData(reader)
            seqid <- match(res@seq_name[1], rd$SN)
            if(!is.na(seqid))
            {
                cat("\rReading file ", format(i, witdh=3))
                # (-1): Correction for 0-based coordinates in alignDepth...
                coords <- c(rd$ID[seqid], (res@coords - 1))
                brg <- bamRange(reader, coords)

                # Add data to global gapSiteList
                ssl <- siteList(brg)
                glsl <- merge(glsl, ssl)

                # Extract gapSiteList from bamRange and merge
                # with existing gapSiteList for appropriate group
                gsl[[groups[i]]] <- merge(gsl[[groups[i]]], ssl)
                bamClose(reader)
                ad <- alignDepth(brg)
                res@ald[, i] <- ad@depth
            }else{
                cat("\n[geneAlignDepth] Missing match for seq_name \"",
                    res@seq_name[1],
                    "\"for file", i, "\n")
            }
        }

        res@gapSites <- as.data.frame(glsl)
        gsl <- lapply(gsl, as.data.frame)
        assign("gsl", gsl, envir=res@ev)

        return(res)
    }
)


setMethod("show", "geneAlignDepth", function(object)
{
    bm<-Sys.localeconv()[7]

    cat("An object of class '", class(object), "'.\n", sep="")
    cat("Length   : ", length(object),      "\n")
    cat("Gene id  : ", object@gene_id[1]  , "\n")
    cat("Gene name: ", object@gene_name[1], "\n")
    cat("Seqid    : ", object@seq_name[1],  "\n")
    cat("Strand   : ", object@strand[1],    "\n")
    cat("nAligns  : ", format(sum(object@nAligns), big.mark=bm), "\n")
    cat("Start    : ", format(object@coords[1],    big.mark=bm), "\n")
    cat("End      : ", format(object@coords[2],    big.mark=bm), "\n")

    if(exists("exons", where=object@ev, inherits=FALSE))
    {
        n<-min(nrow(object@ev$exons), 6L)
        bm<-Sys.localeconv()[7]
        cat("Exon Nr  :",
            format(nrow(object@ev$exons), big.mark = bm),
            "\n")

        print(object@ev$exons[1:n, ])
    }
})

plot.geneAlignDepth <- function(x, col=NULL, ...)
{

    bm<-Sys.localeconv()[7]

    mxy <- max(x@ald)
    length <- ncol(x@ald)
    nGroups <- length(levels(x@group))
    groupNames <- levels(x@group)
    groupv <- as.numeric(x@group)
    nAligns <- tapply(as.numeric(x@nAligns), groupv, sum, na.rm=TRUE)

    # Omits light gray ....
    if(is.null(col))
        col <- terrain.colors(nGroups + 1)[1:nGroups]

    if(length(col)==1)
        col <- rep(col, nGroups)

    if(length(col)!= nGroups)
        stop("One colour for each group (or only one colour) required")

    if(any(nAligns==0))
    {
        warning("0 nAligns present. No align - normalization")
        nAligns <- rep(1, nGroups)
    }else{
        cat("[plot.geneAlignDepth] Renorming to",
            format(nAligns[1], big.mark = bm, width=12),
            "\n"
        )
        nAligns <- nAligns/nAligns[1] # Renorm
    }

    # RPKM: 1e9 * c / ( n * l)
    # c = reads mapped to gene
    # n = total number of reads in sample
    # l = total exon length in gene

    xv <- as.numeric(rownames(x@ald))

    mtx <- matrix(0, nrow=nrow(x@ald), ncol=nGroups)
    rownames(mtx) <- rownames(x@ald) # genomic positions
    colnames(mtx) <- groupNames

    for(i in 1:nGroups)
        mtx[, i] <- apply(as.matrix(x@ald[, x@group==groupNames[i]]), 1, mean)/nAligns[i]
    maxy <- max(mtx)

    plot(xv, mtx[, 1], ylim=c(0,maxy), type="l",
        las=1, bty="n", col=col[1],
        main=x@gene_name,
        xlab=paste("Position on",x@seq_name),
        ylab="Aligns (normalized)", ...)

    if(nGroups > 1)
    {
        for(i in 2:nGroups)
            lines(xv, mtx[, i], col=col[i], ...)
    }

    legend("topright", lwd=2, col=col[1:nGroups],
           legend=groupNames, bty="n")
}


# setGeneric("extractGeneRegions",
#     function(src, trg, gl) standardGeneric("extractGeneRegions"))

setMethod("extractGeneRegions", c("bamReader", "bamWriter", "geneList"),
    function(src, trg, gl)
    {
        n <- length(gl)
        nAligns <- numeric(n)
        for(i in 1:n)
        {
            gm <- gl@l[[i]]
            rng <- readRange(src, gm@seq_name, gm@coords)
            nAligns[i] <- bamSave(trg, rng)
        }
        return(invisible(nAligns))
    }
)


## + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + ##
## CLASS exonAlignDepth
## + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + ##

.exonAlignDepth <- setClass("exonAlignDepth",
    representation(ald="matrix",
        gene_id="character",
        gene_name="character",
        seq_name="character",
        group="factor",
        strand="character",
        label="character",
        nAligns="numeric",
        aldRatio="data.frame",
        junctions="data.frame"
    )
)

setMethod("initialize", "exonAlignDepth", function(.Object)
{
    # Default constructor
    return(.Object)
})

setMethod("length", "exonAlignDepth", function(x){
    return(length(x@group))
})

setMethod("show", "exonAlignDepth", function(object)
{
    bm<-Sys.localeconv()[7]

    cat("An object of class '", class(object), "'.\n", sep="")
    cat("Length   : ", length(object),      "\n")
    cat("Gene id  : ", object@gene_id[1]  , "\n")
    cat("Gene name: ", object@gene_name[1], "\n")
    cat("Seqid    : ", object@seq_name[1],  "\n")
    cat("Strand   : ", object@strand[1],    "\n")
    cat("nAligns  : ", format(sum(object@nAligns), big.mark=bm), "\n")
})



setMethod("[", signature="exonAlignDepth", function(x, i)
{
    i <- sort(unique(as.integer(i)))

    res <- .exonAlignDepth()
    res@ald <- x@ald[,i]
    res@gene_id <- x@gene_id[i]
    res@gene_name <- x@gene_name[i]
    res@seq_name <- x@seq_name[i]
    res@strand <- x@strand[i]
    res@group <- x@group[i]
    res@label <- x@label[i]
    res@nAligns <- x@nAligns[i]
    return(res)
})




setGeneric("exonAlignDepth",
    function(object, ratioLim=5, infVal=1000)
                                    standardGeneric("exonAlignDepth"))


setMethod("exonAlignDepth", "geneAlignDepth",
                    function(object, ratioLim=5, infVal=1000)
{
    res <- .exonAlignDepth()

    # Copy values from incoming object
    res@gene_id <- object@gene_id
    res@gene_name <- object@gene_name
    res@seq_name <- object@seq_name
    res@strand <- object@strand
    res@group <- object@group
    res@label <- object@label
    res@nAligns <- object@nAligns

    # Gene positions
    start <- object@coords[1]
    end   <- object@coords[2]


    ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ##
    ## Identification of Exon junctions
    ## The identification is a two step process:
    ## Putative junction sites are alignment-gap-sites
    ##
    ## A selection criterion for gap-sites an abrupt change in
    ## align-depth (ALD), identified by spikes in ALD ratios
    ## between adjacent positions.
    ## The direction of ALD is selected so that ratios >> 1 are used.
    ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ##

    ## - - - - - - - - - - - - - - - - - - - - ##
    ## A) Using alignment depth ratios
    ## - - - - - - - - - - - - - - - - - - - - ##

    ald <- data.frame(position=as.numeric(rownames(object@ald)),
                      ald=apply(object@ald, 1, mean))

    ## - - - - - - - - - - - - - - - - - - - - ##
    ## Ratio: Look ahead (-> rstart edges)
    ## Relation between in-place and next-position align depth
    ## - - - - - - - - - - - - - - - - - - - - ##
    ald$nxt <- c(ald$ald[-1], ald$ald[nrow(ald)])
    ald$nr  <- ald$ald/ald$nxt

    # Division by zero
    ald$nr[is.infinite(ald$nr)] <- infVal
    # Division 0/0
    ald$nr[is.na(ald$nr)] <- 1

    # Change to integral numbers as we look at ratios >> 1
    ald$nr <- as.integer(ald$nr)

    ## - - - - - - - - - - - - - - - - - - - - ##
    ## Ratio: Go back (-> lend edges)
    ## Relation between in-place and previous-position align depth
    ## - - - - - - - - - - - - - - - - - - - - ##
    ald$prv <- c(ald$ald[1], ald$ald[-nrow(ald)])
    ald$pr <- ald$ald/ald$prv

    # Division by zero
    ald$pr[is.infinite(ald$pr)] <- infVal
    # Division 0/0
    ald$pr[is.na(ald$pr)] <- 1

    # Change to integral numbers as we look at ratios >> 1
    ald$pr <- as.integer(ald$pr)


    ## - - - - - - - - - - - - - - - - - - - - ##
    ## B) Using gap-sites as list of
    ##    putative junctions
    ## - - - - - - - - - - - - - - - - - - - - ##
    junc <- object@gapSites

    mtc <- match(junc$lend, ald$position)
    junc$lend_rat <- ald$nr[mtc]
    mtc <- match(junc$rstart, ald$position)
    junc$rstart_rat <- ald$pr[mtc]

    # Minimum of left and right ratio
    junc$rat <- pmin(junc$lend_rat, junc$rstart_rat)
    junc <- junc[!is.na(junc$rat), ]

    junc$lm_sum <- NULL
    junc$lcl <- NULL
    junc$mcl <- NULL
    junc$lstart <- NULL
    junc$rend <- NULL

    # Copy raw data into returned object
    res@aldRatio <- ald
    res@junctions <- junc


    ## - - - - - - - - - - - - - - - - - - - - ##
    ## C) Set filter and extract
    ##    segmentized ALD matrix
    ## - - - - - - - - - - - - - - - - - - - - ##

    # Filter
    segjunc <- junc[junc$rat > ratioLim, ]

    # Extract segmentated align-depth output matrix
    s <- segmentize(object@ald,
                    begin=segjunc$lend+1,
                    end=segjunc$rstart - 1,
                    start, invert=TRUE)

    # There may be just one column...
    # which then needs to be transformed back into a matrix
    if(is.matrix(s)){
        res@ald <- s
    }else{
        res@ald <- matrix(s)
        rownames(res@ald) <- names(s)
    }

    return(res)
})

setGeneric("aldRatio", function(object) standardGeneric("aldRatio"))
setMethod("aldRatio", "exonAlignDepth", function(object){
    return(object@aldRatio)
})

setGeneric("junctionSites", function(object) standardGeneric("junctionSites"))
setMethod("junctionSites", "exonAlignDepth", function(object){
    return(object@junctions)
})


setGeneric("groupAldMatrix", function(object, renorm=TRUE, f=mean)
                                    standardGeneric("groupAldMatrix"))

setMethod("groupAldMatrix", "exonAlignDepth",
            function(object, renorm=TRUE, f=mean){

    nGroups <- length(levels(object@group))
    groupNames <- levels(object@group)
    groupv <- as.numeric(object@group)
    rawAligns <- tapply(as.numeric(object@nAligns), groupv, sum, na.rm=TRUE)

    if(!is(f, "function"))
        stop("Grouping argument f must be a function")

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
    # Renorming of alignment depth values using absolute read
    # numbers is needed for plotting alignment depth from different
    # samples together
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
    if(any(rawAligns == 0))
    {
        warning("0 nAligns present. No align - normalization")
        nAligns <- rep(1, nGroups)
    }else{
        # Renorm factors
        nAligns <- rawAligns/rawAligns[1]
    }

    mtx <- matrix(0, nrow=nrow(object@ald), ncol=nGroups)
    rownames(mtx) <- rownames(object@ald) # genomic positions
    colnames(mtx) <- groupNames

    for(i in 1:nGroups)
        mtx[, i] <- apply(as.matrix(object@ald[, object@group==groupNames[i]]),
                            1,
                            f) / nAligns[i]


    attr(mtx, "rawAligns") <- rawAligns
    attr(mtx, "nAligns") <- nAligns
    return(invisible(mtx))
})


setGeneric("getNormFactor", function(object) standardGeneric("getNormFactor"))
setMethod("getNormFactor", "exonAlignDepth", function(object){

    groupv <- as.numeric(object@group)
    rawAligns <- tapply(as.numeric(object@nAligns), groupv, sum, na.rm=TRUE)

    if(any(rawAligns == 0))
        return(1)

    return(rawAligns[1])
})

setGeneric("groupAldTable", function(object, renorm=TRUE, f=mean)
                                standardGeneric("groupAldTable"))

setMethod("groupAldTable",  "exonAlignDepth",
            function(object, renorm=TRUE, f=mean)
{

    nGroups <- length(levels(object@group))
    groupNames <- levels(object@group)
    groupv <- as.numeric(object@group)
    rawAligns <- tapply(as.numeric(object@nAligns), groupv, sum, na.rm=TRUE)

    if(!is(f, "function"))
        stop("Grouping argument f must be a function")


    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
    # Renorming of alignment depth values using absolute read
    # numbers is needed for plotting alignment depth from different
    # samples together
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
    if(any(rawAligns == 0))
    {
        warning("0 nAligns present. No align - normalization")
        nAligns <- rep(1, nGroups)
    }else{
        # Renorm factors
        nAligns <- rawAligns/rawAligns[1]
    }

    position <- as.numeric(rownames(object@ald))
    npos <- length(position)

    res <- data.frame(gpos=rep(position, nGroups),
                tpos=1:npos,
                group=rep(groupNames, each=npos),
                ald=numeric(npos * nGroups))

    x <- 1:npos
    for(i in 1:nGroups)
    {
        res$ald[x] <- apply(as.matrix(object@ald[,
                                object@group==groupNames[i]]),
                            1, f) / nAligns[i]

        x <- x + npos
    }

    attr(res, "rawAligns") <- rawAligns
    attr(res, "nAligns") <- nAligns

    return(invisible(res))
})


## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ##
## Cut out regions with low alignment coverage
## e.g. due to previously undetected introns
## The function is not restricted to contiguous areas!
## Therefore junction positions are omitted
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ##

setGeneric("cutFlatAlignDepth", function(object, ratio=0.02, f=mean)
    standardGeneric("cutFlatAlignDepth"))

setMethod("cutFlatAlignDepth", "exonAlignDepth",
            function(object, ratio=0.02, f=mean)
{

    if(!is(f, "function"))
        stop("Grouping argument f must be a function")

    ratio <- ratio[1]
    if(ratio <= 0)
        stop("ratio must be >= 0!")

    if(ratio >= 0.5)
        stop("ratio must be < 0.5!")

    res <- object
    mam <- groupAldMatrix(object, f=f)
    alm <- apply(mam, 1, mean)
    lim <- max(alm, na.rm=TRUE) * ratio
    alm[alm < lim] <- 0

    # Cut out all regions below limit
    m <- res@ald[alm > 0, ]

    # object@ald may only have one column ...
    if(!is.matrix(m)){
        m <- matrix(m, ncol=ncol(object@ald))
        rownames(m) <- rownames(object@ald)[alm > 0]
    }
    res@ald <- m

    # Empty junctions table
    res@junctions <- res@junctions[0]

    return(res)
})



plot.exonAlignDepth <- function(x,
        col=NULL, ylim=NULL, xlim=NULL, xunit=1000, yunit=1000, ...)
{

    bm<-Sys.localeconv()[7]

    nGroups <- length(levels(x@group))
    groupNames <- levels(x@group)
    groupv <- as.numeric(x@group)
    nAligns <- tapply(as.numeric(x@nAligns), groupv, sum, na.rm=TRUE)

    if(is.null(col))
    {
        # Omits light gray
        col <- terrain.colors(nGroups + 1)[1:nGroups]
    }else if(length(col)==1){
        col <- rep(col, nGroups)
    }else if(length(col)!= nGroups){
        stop("One colour for each group (or only one colour) required")
    }


    if(any(nAligns==0))
    {
        warning("0 nAligns present. No align - normalization")
        nAligns <- rep(1, nGroups)
    }else{
        cat("[plot.geneAlignDepth] Renorming to",
            format(nAligns[1], big.mark = bm, width=12),
            "\n"
        )
        nAligns <- nAligns/nAligns[1] # Renorm
    }

    xlen <- nrow(x@ald)
    xv <- 1:xlen

    mtx <- matrix(0, ncol=nrow(x@ald), nrow=nGroups)
    colnames(mtx) <- rownames(x@ald) # genomic positions
    rownames(mtx) <- groupNames

    for(i in 1:nGroups)
        mtx[i, ] <- apply(as.matrix(x@ald[, x@group==groupNames[i]]),
                                        1,
                                        mean) / nAligns[i]

    mxy <- max(mtx)

    # Units for scaling x and y axis
    xunit <- '^'(10, max(floor(log10(xlen)) - 1, 0))
    yunit <- '^'(10, max(floor(log10(mxy)) - 1, 0))


    # Scale x and y axis
    if(is.null(xlim))
        xlim <- c(0, ceiling(xlen/xunit) * xunit)

    if(is.null(ylim))
        ylim <- c(0, ceiling(mxy/yunit) * yunit)

    plot(xv, mtx[1, ],
        ylim=ylim,
        xlim=xlim,
        type="l",
        las=1,
        bty="n", col=col[1],
        main=x@gene_name,
        xlab=paste("Position on",x@seq_name),
        ylab="Aligns (normalized)", ...)

    if(nGroups > 1)
    {
        for(i in 2:nGroups)
            lines(xv, mtx[i, ], col=col[i], ...)
    }

    legend("topright", lwd=2, col=col[1:nGroups],
           legend=groupNames, bty="n")
}


## + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + ##
##
## + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + ##

.exonLoessModel <- setClass("exonLoessModel",
        representation(gene_id="character",
            gene_name="character",
            seq_name="character",
            group="character",
            strand="character",
            nAligns="numeric",
            loessPred="matrix" # Loess predicted values
        )
)

setMethod("initialize", "exonLoessModel", function(.Object)
{
    # Default constructor
    return(.Object)
})

setMethod("length", "exonLoessModel", function(x){
    return(length(x@group))
})

setMethod("show", "exonLoessModel", function(object)
{
    bm<-Sys.localeconv()[7]

    cat("An object of class '", class(object), "'.\n", sep="")
    cat("Length   : ", length(object),      "\n")
    cat("Gene id  : ", object@gene_id[1]  , "\n")
    cat("Gene name: ", object@gene_name[1], "\n")
    cat("Seqid    : ", object@seq_name[1],  "\n")
    cat("Strand   : ", object@strand[1],    "\n")
    cat("nAligns  : ", format(sum(object@nAligns), big.mark=bm), "\n")
})



setMethod("[", signature="exonLoessModel", function(x, i)
{
    i <- sort(unique(as.integer(i)))

    res <- .exonAlignDepth()
    res@loessPred <- x@loessPred[,i]
    res@gene_id <- x@gene_id[i]
    res@gene_name <- x@gene_name[i]
    res@seq_name <- x@seq_name[i]
    res@strand <- x@strand[i]
    res@group <- x@group[i]
    res@nAligns <- x@nAligns[i]
    return(res)
})




setGeneric("exonLoessModel", function(object, span=0.1, f=mean)
                                standardGeneric("exonLoessModel"))

setMethod("exonLoessModel", "exonAlignDepth",
    function(object, span=0.1, f=mean)
    {
        ald <- object@ald
        npos <- nrow(ald)

        # Sample groups
        group <- object@group
        grl <- levels(group)
        ngrp <- length(grl)

        res <- .exonLoessModel()
        res@gene_id <- object@gene_id
        res@gene_name <- object@gene_name
        res@seq_name <- object@seq_name
        res@strand <- object@strand
        res@nAligns <- object@nAligns

        res@loessPred <-  matrix(0, nrow=npos, ncol=ngrp)
        rownames(res@loessPred) <- rownames(object@ald)

        # Change group assignment: One column per group
        mam <- groupAldMatrix(object, f=f)
        res@group <- colnames(mam)
        colnames(res@loessPred) <- colnames(mam)

        x <- 1:nrow(mam)
        for(i in 1:ncol(mam))
        {
            lmod <- loess(mam[, i]~x, span=span)
            res@loessPred[, i] <- predict(lmod, x, se=FALSE)
        }
        return(res)
})

setMethod("cutFlatAlignDepth", "exonLoessModel", function(object, ratio=0.02){

    ratio <- ratio[1]
    if(ratio <= 0)
        stop("ratio must be >= 0!")

    if(ratio >= 0.5)
        stop("ratio must be < 0.5!")

    res <- object
    mam <- object@loessPred

    # Loess model is flat line
    if(max(mam==0))
        return(res)

    alm <- apply(mam, 1, mean)
    lim <- max(alm, na.rm=TRUE) * ratio
    alm[alm < lim] <- 0

    # Cut out all regions below limit
    m <- res@loessPred[alm > 0, ]

    # object@loessPred may only have one column ...
    if(!is.matrix(m)){
        m <- matrix(m, ncol=ncol(object@loessPred))
        rownames(m) <- rownames(object@loessPred)[alm > 0]
    }
    res@loessPred <- m

    return(res)
})



plot.exonLoessModel <- function(x, col=NULL, lwd=2, ...)
{
    ng <- length(x)

    # Omits light gray ....
    if(is.null(col))
        col <- terrain.colors(ng + 1)[1:ng]

    if(length(col)==1)
        col <- rep(col, ng)

    if(length(col) != ng)
        stop("One colour per group required!")

    lspr <- x@loessPred
    # Negative align depth seems to be implausible
    lspr[lspr < 0] <- 0

    mxy <- max(lspr)
    xlen <- nrow(lspr)
    xv <- 1:xlen

    # Units for scaling x and y axis
    xunit <- '^'(10, max(floor(log10(xlen)) - 1, 0))
    yunit <- '^'(10, max(floor(log10(mxy)) - 1, 0))

    # Scale axes using scale units
    yul <- ceiling(mxy/yunit) * yunit
    xul <- ceiling(nrow(lspr)/xunit) * xunit

    plot(xv, lspr[, 1],
            xlim=c(0, xul), ylim=c(0, yul),
            las=1, bty="n", type="l",
            col=col[1], lwd=lwd,
            xlab= paste("Position (seq=", x@seq_name, ")", sep=""),
            ylab="Alignment Depth")

    if(ng > 1)
    {
        for(i in 2:ng)
            lines(xv, lspr[, i], col=col[i], lwd=lwd)
    }

    legend("topright", legend=x@group, col=col[1:3], lwd=lwd, bty="n")

    title(main=x@gene_name,
          sub=paste("Gene id: ", x@gene_id,
                    ", (", x@strand, ") - Strand", sep=""))
}


setGeneric("groupRatio",
    function(object, lim=1.2, cut=0, order=NULL, f=mean)
                                    standardGeneric("groupRatio"))


setMethod("groupRatio", "exonLoessModel",
        function(object, lim=1.2, cut=0, order=NULL, f=mean)
{

    if(!is(f, "function"))
        stop("Grouping argument f must be a function")

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
    # We possibly remove areas with low align depth
    # e.g. long regions due to introns
    # in order to increase sensitivity
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
    cut <- cut[1]
    if(cut > 0)
        object <- cutFlatAlignDepth(object, ratio=cut)

    if(!is.null(order))
    {
        if(!all(is.element(order, 1:ncol(object@loessPred))))
            stop("order values out of range!")
    }else{
        order <- 1:ncol(object@loessPred)
    }

    lim <- lim[1]
    if(lim <= 1)
        stop("lim must be greater than 1")

    ng <- length(object)
    if(ng==1)
        return(1)

    lspr <- object@loessPred
    lspr[lspr < 0] <- 0

    lspr <- lspr[, order]

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
    # Two ratios are calculated so that ascending and descending ratios
    # can be identified using the same ratio-limit
    # ascrat: Ascending ratio
    # decrat: Descending ratio
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
    ascrat <- matrix(0, nrow=nrow(lspr), ncol=ng - 1)
    decrat <- matrix(0, nrow=nrow(lspr), ncol=ng - 1)
    for(i in 1:(ng - 1))
    {
        ascrat[, i] <- lspr[, i + 1] / lspr[, i]
        decrat[, i] <- lspr[, i] / lspr[, i + 1]
    }

    # Remove unwanted extreme values
    ascrat[is.na(ascrat)] <- 0
    ascrat[is.nan(ascrat)] <- 0
    ascrat[is.infinite(ascrat)] <- 0

    decrat[is.na(decrat)] <- 0
    decrat[is.nan(decrat)] <- 0
    decrat[is.infinite(decrat)] <- 0


    ascmin <- apply(ascrat, 1, min)
    asc_lim <- sum(ascmin > lim) / length(ascmin)

    decmin <- apply(decrat, 1, min)
    dec_lim <- sum(decmin > lim) / length(decmin)

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
    # Returned value will always be in ascending view:
    # Ascending ratio: >1
    # Descening ratio: <1
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
    return(ifelse(asc_lim > dec_lim, asc_lim, (-1) * dec_lim))
})


## + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + ##
##
## + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + ##

setGeneric("saveAldData",
    function(bs, gl, path, order=NULL, lim=1.1, startId=1, f=mean)
        standardGeneric("saveAldData"))

setMethod("saveAldData", c("sampleBamFiles", "geneList"),
    function(bs, gl, path, order=NULL, lim=1.1, startId=1, f=mean)
{
    mb <- "[saveAldData]"

    groups <- sampleGroups(bs)
    nGroups <- length(levels(groups))
    if(!is.null(order))
    {
        if(!all(is.element(order, 1:nGroups)))
            stop("order values out of range!")
    }else{
        order <- 1:nGroups
    }

    gadp <- file.path(path, "gene_align_depth")
    eadp <- file.path(path, "exon_align_depth")
    elom <- file.path(path, "exon_loess_model")
    celom <- file.path(path, "cut_exon_loess_model")

    if(!file.exists(path))
        dir.create(path)

    if(!file.exists(gadp))
        dir.create(gadp)
    if(!file.exists(eadp))
        dir.create(eadp)
    if(!file.exists(elom))
        dir.create(elom)
    if(!file.exists(celom))
        dir.create(celom)

    cat(mb, "nGenes: ", length(gl), "\n")

    n <- length(gl)

    aldrat <- data.frame(id=1:n, gene_name="", maxald=0, gr=0, cgr=0,
                      filename="", stringsAsFactors=FALSE)

    for(i in startId:n)
    {
        gm <- gl[i]

        cat(mb, "i: ", format(i, width=3), "\t gene: ", gm@gene_name, "\n")

        gad <- geneAlignDepth(bs, gm)
        ead <- exonAlignDepth(gad, ratioLim=5)
        elm <- exonLoessModel(ead, span=0.1, f=f)
        celm <- cutFlatAlignDepth(elm, ratio=0.1, f=f)

        aldrat$maxald[i] <- max(gad@ald)

        aldrat$gene_name[i] <- gm@gene_name
        aldrat$gr[i] <- groupRatio(elm, order=order, lim=lim, f=f)
        aldrat$cgr[i] <- groupRatio(celm, order=order, lim=lim, f=f)

        cat("\n")

        filename <- paste(i, "_", gm@gene_name, ".pdf", sep="")
        aldrat$filename[i] <- filename

        pdf(file.path(gadp, filename))
        plot(gad)
        dev.off()

        pdf(file.path(eadp, filename))
        plot(ead)
        dev.off()

        pdf(file.path(elom, filename))
        plot(elm)
        dev.off()

        pdf(file.path(celom, filename))
        plot(celm)
        dev.off()

        write.table(aldrat, file=file.path(path, "aldratio.csv"),
                    sep=";", dec=",", row.names=FALSE)
    }
    save(aldrat, file=file.path(path, "aldratio.RData"))
    cat(mb, "Finished.\n")
    return(invisible(aldrat))
})


## + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + ##
## extractGeneRegions:
## Copy alignments from genetic regions out to second BAM files
## for multiple files at once (using sampleBamFiles)
## + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + ##

setMethod("extractGeneRegions", c("bamReader", "character", "geneList"),
    function(src, trg, gl)
    {
        writer <- bamWriter(getHeader(src), trg[1])
        nAligns <- sum(extractGeneRegions(src, writer, gl))
        bamClose(writer)
        return(invisible(nAligns))
    }
)

setMethod("extractGeneRegions", c("sampleBamFiles", "sampleBamFiles", "geneList"),
    function(src, trg, gl)
    {
        n <- length(src)
        if(length(trg) < n)
            stop("Not enough targets for writing (length(sampleBamFiles))")

        cat("[extractGeneRegions] Writing", n, "files.\n")

        nAligns <- numeric(n)
        for(i in 1:n)
        {
            cat("[extractGeneRegions] File ", i, ".\n")
            reader <- bamReader(bamFiles(src)[i], bamIdxFiles(src)[i])
            nAligns[i] <- extractGeneRegions(reader, bamFiles(trg)[i], gl)
            bamClose(reader)
        }
        return(invisible(nAligns))
    }
)


setMethod(f="bamSort",
    signature=c("sampleBamFiles"),
    definition=function(object,
        prefix="sorted",
        byName=FALSE,
        maxmem=1e+9,
        path=dirname(bamFiles(object)))
    {
        byName <- byName[1]
        maxmem <- floor(maxmem[1])
        n <- length(object)

        if(length(prefix)==1 & n > 1)
            prefix <- paste(prefix, 1:n, sep="_")

        if(length(prefix)!=n)
            stop("There must be one prefix for each file (or one prefix)")

        message("[bamSort] Maxmem  : ", maxmem)
        message("[bamSort] By Name : ", byName)


        for(i in 1:n)
        {
            cat("\r[bamSort] i:", format(i, w=3) , "/", n, ".")
            reader <- bamReader(bamFiles(object)[i])

            .Call(C_bam_reader_sort_file,
                reader@filename,
                path.expand(file.path(path, prefix[i])),
                maxmem, byName, PACKAGE="rbamtools")

        }
        cat("\n[bamSort] Sorting finished.\n")
        return(invisible(paste(prefix, "bam", sep=".")))
    }
)

## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ##
##   Index related functions
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ##

setMethod("createIndex",
    c("sampleBamFiles", "character"),
    function(object, idx_filename)
    {
        n <- length(object)
        if(length(idx_filename)!=n)
            stop("There must be one index file name for each BAM file!")

        cat("[createIndex] sampleBamFiles n: ", n, "\n")
        for(i in 1:n)
        {
            cat("\r[createIndex] ", format(i, width=3))
            .Call(C_bam_reader_create_index,
                path.expand(bamFiles(object)[i]),
                path.expand(idx_filename[i]), PACKAGE="rbamtools")
        }
        cat("\n[createIndex] Finished.\n")
    }
)

setMethod("createIndex",
    c("sampleBamFiles", "missing"),
    function(object, idx_filename)
    {
        n <- length(object)

        cat("[createIndex] sampleBamFiles n: ", n, "\n")
        for(i in 1:n)
        {
            cat("\r[createIndex] ", format(i, width=3))
            .Call(C_bam_reader_create_index,
                  path.expand(bamFiles(object)[i]),
                  path.expand(bamIdxFiles(object)[i]), PACKAGE="rbamtools")
        }
        cat("\n[createIndex] Finished.\n")
    }
)


## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ##
## END OF FILE
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ##
