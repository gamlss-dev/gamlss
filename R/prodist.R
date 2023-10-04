## S3 method for extracting fitted/predicted distributions3 objects
## associated methods are in gamlss.dist (as well as distributions3, topmodels, etc.)
prodist.gamlss <- function(object, ...) {
  d <- predictAll(object, ...)
  d$y <- NULL
  class(d) <- c("GAMLSS", "distribution")
  return(d)  
}
