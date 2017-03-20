

multSeq <- function(beg, end)
{
    return(.Call(C_mult_seq,
            as.integer(beg), as.integer(end), PACKAGE="rbamtools"))
}




## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ##
## Segmentation
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ##

setGeneric("segmentize",
    function(x, begin, end,
        offset=1, margin=1, invert=FALSE) standardGeneric("segmentize"))


setMethod("segmentize", "ANY",
    function(x, begin, end, offset=1, margin=1, invert=FALSE)
    {
        if(any(end < begin))
            stop("end < begin found!")

        if(length(begin) != length(end))
            stop("begin and end must have equal length!")

        sb <- as.integer(begin - offset[1] + 1)
        se <- as.integer(end - offset[1] + 1)


        index <- multSeq(sb, se)

        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
        # Eventually switch to complement of selected segments
        # Then, ordering cannot be changed...
        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
        if(invert[1])
        {
            lg <- rep(TRUE, length(x))
            lg[index] <- FALSE
            index <- lg
        }

        return(x[index])
    }
)


setMethod("segmentize", "matrix",
    function(x, begin, end, offset=1, margin=1, invert=FALSE)
    {
        if(any(end < begin))
            stop("end < begin found!")

        if(length(begin) != length(end))
            stop("begin and end must have equal length!")

        sb <- as.integer(begin - offset[1] + 1)
        se <- as.integer(end - offset[1] + 1)

        index <- multSeq(sb, se)

        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
        # Eventually switch to complement of selected segments
        # Then, ordering cannot be changed...
        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #

        if(margin[1] == 1)
            len <- nrow(x)
        else
            len <- ncol(x)

        if(invert[1])
        {
            lg <- rep(TRUE, len)
            lg[index] <- FALSE
            index <- lg
        }


        if(margin[1] == 1)
            return(x[index, ])
        else
            return(x[, index])
    }
)

setMethod("segmentize", "data.frame",
    function(x, begin, end, offset=1, margin=1, invert=FALSE)
    {
        if(any(end < begin))
            stop("end < begin found!")

        if(length(begin) != length(end))
            stop("begin and end must have equal length!")

        sb <- as.integer(begin - offset[1] + 1)
        se <- as.integer(end - offset[1] + 1)

        index <- multSeq(sb, se)

        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
        # Eventually switch to complement of selected segments
        # Then, ordering cannot be changed...
        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #

        if(margin[1] == 1)
            len <- nrow(x)
        else
            len <- ncol(x)

        if(invert[1])
        {
            lg <- rep(TRUE, len)
            lg[index] <- FALSE
            index <- lg
        }


        if(margin[1] == 1)
            return(x[index, ])
        else
            return(x[, index])
    }
)




## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ##
## END OF FILE
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ##
