
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ##
## Test vector segmentation
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ##

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

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# END OF FILE
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
