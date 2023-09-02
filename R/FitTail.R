# library(gamlss)
# those are function for checking the tails of a distribution
# created by Mikis Stasinopoulos Bob Rigby and Vlasios Voudouris Dec 2011 
# rewriten MARCH 2018 by Mikis
# TO DO
# 
################################################################################
################################################################################
################################################################################
################################################################################
# log log Survival plots
# there are three functions for that plus one which connect them al
# i)   loglogSurv1()  : method 2 for for type I tails (log WEI)
# ii)  loglogSurv2()  : method 2 for type II tails (WEI)
# iii) loglogSurv3()  : method 2 for type III tails (GUMBEL)
# iv)  loglogSurv     : combines the above functions and selects the best
################################################################################
# log Survival plot
# v)    logSurv() : plots the empirical survival function log(1-ecdf) 
#        against log(y) for part or all of the data and fits three linear model
# vi)  logSurvo()  plots the empirical survival function log(1-ecdf) or ecdf 
#        against log(y) for part or all the data (no fitting here)
# vii)  loglogplot() plots the empirical survival function log(1-ecdf) for all data
# viii)  ECDF()
#       
#-------------------------------------------------------------------------------
# those are in gamlss.tr
# vii)  fitTail : fits a truncated gamlss.family distribution to the tail of the data
# viii) fitTailAll : fits a (Hill type) series of Fit using the fitTail function 
################################################################################
################################################################################
################################################################################
################################################################################
# TYPE I
################################################################################
loglogSurv1 <- function(y, 
                         prob = 0.9, 
                         print = TRUE,
                         title = NULL,
                          lcol = gray(.1),
                         ltype = 1,
                        ...)
{
#-----------------
  Xlab1 <- paste("log(log(", deparse(substitute(y)), "))", sep="")
    cdF <- ECDF(y) # get ecdf for all data
     mY <- quantile(y, probs=prob)
      Y <- y[y>mY]
if (mY<1) stop("For the method for Type I to work the minimum value of the tail must be greater than 1") 
   newY <- log(-log(1-cdF(Y)))
  place <- "bottomleft"
   Ylab <- paste("log(1-F(",deparse(substitute(y)), "))", sep="")
howmany <- length(Y)
      x <- log(log(Y))
     m1 <- gamlss(newY ~ x, trace=FALSE)
    ess <- sum((newY-fitted(m1))^2)
   Ylab <- paste("log(-log( S(", deparse(substitute(y)), ")))", sep="")  
  if (print)
  {
    cat("coefficients",  coef(m1), "\n")
    cat("error sum of squares", ess, "\n")
  }
  plot(newY~x, xlab=Xlab1, ylab=Ylab, ...)
  lines(fitted(m1)~x, col=lcol, lty=ltype)  
  if (is.null(title))
    {
    title(paste("Log Log Survival plot (Type I) for",  (1-prob)*100,
                          "% of the right tail of",  deparse(substitute(y))))
    } else title(title)
  invisible(m1)  
}
################################################################################
################################################################################
################################################################################
################################################################################
# Type II
################################################################################
loglogSurv2 <- function(y, 
                          prob = 0.9, 
                         print = TRUE,
                         title = NULL,
                          lcol = gray(.1),
                         ltype = 1,
                         ...)
{
  #-----------------
  Xlab1 <- paste("log(", deparse(substitute(y)), ")", sep="")
    cdF <- ECDF(y) # get ecdf for all data
     mY <- quantile(y, probs=prob)
      Y <- y[y>mY]
  if (mY<1) stop("For the method for Type I to work the minimum value of the tail must be greater than 1") 
   newY <- log(-log(1-cdF(Y)))
  place <- "bottomleft"
   Ylab <- paste("log(1-F(",deparse(substitute(y)), "))", sep="")
howmany <- length(Y)
      x <- log(Y)
     m1 <- gamlss(newY ~ x, trace=FALSE)
    ess <- sum((newY-fitted(m1))^2)
   Ylab <- paste("log(-log( S(", deparse(substitute(y)), ")))", sep="")  
  if (print)
  {
    cat("coefficients",  coef(m1), "\n")
    cat("error sum of squares", ess, "\n")
  }
  plot(newY~x, xlab=Xlab1, ylab=Ylab, ...)
  lines(fitted(m1)~x, col=lcol, lty=ltype)  
  if (is.null(title))
  {
    title(paste("Log Log Survival plot (Type II) for", (1-prob)*100,
                "% of the right tail of",  deparse(substitute(y))))
  } else title(title)
  invisible(m1)  
}
################################################################################
################################################################################
################################################################################
################################################################################
# Type III
################################################################################
loglogSurv3 <- function(y, 
                     prob = 0.9, 
                    print = TRUE,
                    title = NULL,
                     lcol = gray(.1),
                    ltype = 1,
                        ...)
{
  #-----------------
    Xlab1 <- paste( deparse(substitute(y)), sep="")
      cdF <- ECDF(y) # get ecdf for all data
       mY <- quantile(y, probs=prob)
        Y <- y[y>mY]
  if (mY<1) stop("For the method for Type I to work the minimum value of the tail must be greater than 1") 
     newY <- log(-log(1-cdF(Y)))
    place <- "bottomleft"
     Ylab <- paste("log(1-F(",deparse(substitute(y)), "))", sep="")
  howmany <- length(Y)
        x <- Y
       m1 <- gamlss(newY ~ x, trace=FALSE)
      ess <- sum((newY-fitted(m1))^2)
     Ylab <- paste("log(-log( S(", deparse(substitute(y)), ")))", sep="")  
  if (print)
  {
    cat("coefficients",  coef(m1), "\n")
    cat("error sum of squares", ess, "\n")
  }
  plot(newY~x, xlab=Xlab1, ylab=Ylab, ...)
  lines(fitted(m1)~x, col=lcol, lty=ltype)  
  if (is.null(title))
  {
    title(paste("Log Log Survival plot (Type III) for",  (1-prob)*100, 
                "% of the right tail of",  deparse(substitute(y))))
  } else title(title)
  invisible(m1)  
}
################################################################################
################################################################################
################################################################################
################################################################################
# Select the best from type I II and III
################################################################################
loglogSurv <- function(y, 
                 prob = 0.9, 
                print = TRUE,
                title = NULL,
                 lcol = gray(.1),
                ltype = 1,
                 plot = TRUE,
                      ...)
{
#-----------------
# body of the function starts here
#  Xlab1 <- paste("log(", deparse(substitute(y)), ")", sep="")
  cdF <- ECDF(y) # get ecdf for all data
   mY <- quantile(y, probs=prob)
    Y <- y[y>mY]
  if (mY<1) stop("For the method for Type I to work the minimum value of the tail must be greater than 1") 
 newY <- log(-log(1-cdF(Y)))
 Ylab <- paste("log(1-F(",deparse(substitute(y)), "))", sep="")
  howmany <- length(Y)
    x1 <- log(log(Y))
    m1 <- gamlss(newY ~ x1, trace=FALSE)
  ess1 <- sum((newY-fitted(m1))^2)
    x2 <- log(Y)
    m2 <- gamlss(newY ~ x2, trace=FALSE)
  ess2 <- sum((newY-fitted(m2))^2)
    x3 <- Y
    m3 <- gamlss(newY ~ x3, trace=FALSE)
  ess3 <- sum((newY-fitted(m3))^2) 
   ess <- c(ess1, ess2, ess3)
   num <- which.min(ess)
  matcoef <- rbind(coef(m1), coef(m2), coef(m3))
  matcoef <- cbind(matcoef, ess)
 dimnames(matcoef) <- list(c("type I", "type II", "type III"), c(" Intercept",
" slope", " Error SS"))
if  (print)
{
  cat("Linear regression coefficients", "\n")
  printCoefmat(matcoef, digits = 6, signif.stars = TRUE)
  matk <- matcoef[,-3]
  matk[,1] <-exp(matk[,1])
  colnames(matk) <- c("k:2,4,6", "k:1,3,5")
  cat("Estimates for parameters k", "\n")
  printCoefmat(matk, digits = 3, signif.stars = TRUE) 
}
if (plot)
 {
  
  Ylab <- paste("log(-log( S(", deparse(substitute(y)), ")))", sep="")
  switch(num, 
         { Xlab1 <- paste("log(log(", deparse(substitute(y)), "))", sep="")
           plot(newY~x1, xlab=Xlab1, ylab=Ylab, ...)
           lines(fitted(m1)~x1, col=lcol, lty=ltype)
           if (is.null(title))
           {
             title(paste("Log Log Survival plot (Type I) for", (1-prob)*100, 
                         "% of the right tail of",  deparse(substitute(y))))
           } else title(title)
           },
         {
           Xlab1 <- paste("log(", deparse(substitute(y)), ")", sep="")
           plot(newY~x2, xlab=Xlab1, ylab=Ylab, ...)
           lines(fitted(m2)~x2, col=lcol, lty=ltype)
           if (is.null(title))
           {
             title(paste("Log Log Survival plot (Type II) for", (1-prob)*100, 
                         "% of the right tail of",  deparse(substitute(y))))
           } else title(title)
         },
         {
           Xlab1 <- paste(deparse(substitute(y)), sep="")
           plot(newY~x3, xlab=Xlab1, ylab=Ylab,...)
           lines(fitted(m3)~x3, col=lcol, lty=ltype)
           if (is.null(title))
           {
             title(paste("Log Log Survival plot (Type III) for", (1-prob)*100, 
                         "% of the right tail of",  deparse(substitute(y))))  
           } else title(title)   
           
         }
         )   
 }
model <- switch(num, m1,m2,m3)   
return(model) 
}
################################################################################
################################################################################
################################################################################
################################################################################
# this plots the empirical log(1-ecdf) agains log(y)
# The complementary cumulative distribution function (CCDF) plot
# fits also a linear and a quadraitic fit to the resulting plot yo help 
# interpretation 
################################################################################
logSurv <- function(y, 
                  prob = 0.9, 
                  tail = c("right", "left"), 
                  plot = TRUE, 
                 lines = TRUE,
                 print = TRUE,
                 title = NULL,
                  lcol = c(gray(.1),gray(.2), gray(.3)), 
                 ltype = c(1,2,3),
...)
{
# body of the function starts here
#  require(gamlss)
  tail <- match.arg(tail)
  Xlab <- paste("log(", deparse(substitute(y)), ")", sep="")
   cdF <- ECDF(y) # get ecdf for all data
    mY <- quantile(y, probs=prob)
  if (tail=="right")
        {  Y <- y[y>mY]
        newY <-  log(1-cdF(Y))
       newY2 <- log(-log(1-cdF(Y)))
       place <- "bottomleft"
        Ylab <- paste("log(1-F(",deparse(substitute(y)), "))", sep="")
     howmany <- length(Y)
        } else
        {
           Y <- y[y<mY]
        newY <-  log(cdF(Y))
       newY2 <- log(-log(cdF(Y)))
       place <- "topleft"
        Ylab <- paste("log(F(",deparse(substitute(y)), "))", sep="")
     howmany <- length(Y)
        }
# model fitting
if (lines)
{  #offset(rep(min(newY), length(newY)))
     wm <- which.min(log(Y))
   varY <- newY-newY[wm]
   varX <- log(Y)-log(Y)[wm]
  varX2 <- (log(Y))^2 -(log(Y)[wm])^2
     m1 <- gamlss(varY ~ varX-1, trace=FALSE) # fitting k1=1 model
  #points(fitted(m1)~varX, col="red")
    m2 <- gamlss(varY ~ varX + varX2 - 1, trace=FALSE)  # fitting k1=2 model
  #points(fitted(m2)~varX, col="red")
    m3 <- gamlss(newY2 ~ log(Y), trace=FALSE) # fitting k
  fv3  <- -exp(fitted(m3)) 
}
  if (plot)
  {
    plot(newY[order(Y)]~log(Y)[order(Y)], xlab=Xlab, ylab=Ylab, ...)
    if (lines){
    lines(I(fitted(m1)[order(Y)]+newY[wm])~log(Y)[order(Y)], col=lcol[1], lty=ltype[1], lwd=2)
    lines(I(fitted(m2)[order(Y)]+newY[wm])~log(Y)[order(Y)], col=lcol[2], lty=ltype[2] , lwd=2)
    lines(fv3[order(Y)]~log(Y)[order(Y)], col=lcol[3], lty=ltype[3], lwd=2 )
    }
    if (is.null(title))
    {
      if (tail=="right")
      {
        title(main=paste(paste((1-prob)*100,"%",sep=""),  tail, "tail", "of", 
                         deparse(substitute(y)), ",",howmany,"obs."))
      } else
      {
        title(main=paste(paste((prob)*100,"%",sep=""),  tail, "tail", "of", 
                         deparse(substitute(y)), ",",howmany,"obs."))
      }  
    } else title(title)
    if (lines)
    {
      legend(place, legend=c("linear", "quadratic", "exponential"), col=lcol, lty=ltype, lwd=2 )
    }
  }
}
################################################################################
################################################################################
################################################################################
################################################################################
# THE SAME AS logSurv() WITHOUT FITTING MODELS 
logSurv0 <- function(y, 
                     prob = 0.9, 
                     tail = c("right", "left"), 
                     plot = TRUE, 
                    title = NULL,
                     ...)
{ # body of the function starts here
  #  require(gamlss)
   tail <- match.arg(tail)
   Xlab <- paste("log(", deparse(substitute(y)), ")", sep="")
    cdF <- ECDF(y) # get ecdf for all data
     mY <- quantile(y, probs=prob)
  if (tail=="right")
  {   Y <- y[y>mY]
   newY <-  log(1-cdF(Y))
  newY2 <- log(-log(1-cdF(Y)))
  place <- "bottomleft"
   Ylab <- paste("log(1-F(",deparse(substitute(y)), "))", sep="")
howmany <- length(Y)
  } else
  {
      Y <- y[y<mY]
   newY <-  log(cdF(Y))
  newY2 <- log(-log(cdF(Y)))
  place <- "topleft"
   Ylab <- paste("log(F(",deparse(substitute(y)), "))", sep="")
howmany <- length(Y)
  }
  if (plot)
  {
    plot(newY[order(Y)]~log(Y)[order(Y)], xlab=Xlab, ylab=Ylab, ...)
    if (is.null(title))
    {
      if (tail=="right")
      {
        title(main=paste(paste((1-prob)*100,"%",sep=""),  tail, "tail", "of", 
                         deparse(substitute(y)), ",",howmany,"obs."))
      } else
      {
        title(main=paste(paste((prob)*100,"%",sep=""),  tail, "tail", "of", 
                         deparse(substitute(y)), ",",howmany,"obs."))
      }  
     
    } else title(title)
  }
  M <-rbind(newY[order(Y)], log(Y)[order(Y)])
  rownames(M) <- c("logS", "logy")
  invisible(M) 
}
################################################################################
################################################################################
################################################################################
################################################################################
loglogplot <- function(y, nplus1=TRUE, ...)
{
  Xlab <- paste("log(", deparse(substitute(y)), ")", sep="")
  Ylab <- paste("log(1-F(",deparse(substitute(y)), "))", sep="")   
  if (any(y<=0)) {
     y <- y+abs(min(y))+1
    warning("negative values in y, it is shifted to y+abs(min(y))+1 ")
  }
     Y <- unique(sort(y))
     n <- length(y)
  ecdf <- if (nplus1)  cumsum(tabulate(match(y, Y)))/(n+1)
          else         cumsum(tabulate(match(y, Y)))/n
     y <- if (nplus1) 1-ecdf else  1-c(0, ecdf[-n])
  plot(Y, y, log="xy", xlab=Xlab,ylab=Ylab,...)
  invisible(data.frame(x=Y,y=y))
} 
#------------------------------------------------------------
loglogplot0 <- function(x,  ...)
{
  ff <- ecdf(x)
   x <- unique(sort(x))
   F <- ff(x)
  FF <- c(0, F[-length(F)])
   y <- 1-FF
 plot(x, y, log="xy",...)
   M <- data.frame(x,y)
  invisible(M)
} 
################################################################################
################################################################################
################################################################################
################################################################################
# a function for ecdf the the difference that it divide by n+1
# Mikis 9-1-2010 added weights
################################################################################
################################################################################
ECDF <- function (y, weights=NULL)
{
if (is.null(weights))
  {
  ysort <- unique(sort(y))
      n <- length(y)
   ecdf <- cumsum(tabulate(match(y, ysort)))/(n + 1)
    fun <- stepfun(ysort, c(0, ecdf))
  class(fun) <- c("ecdf", "stepfun")
  return(fun)
  } else
  {
    if (length(weights)!=length(y)) stop("y and weights have different length")
    i <- is.na(weights) | weights == 0
    if (any(i))
    {
      y <- y[!i]
weights <- weights[!i]
    }
  ysort <- unique(sort(y))
      n <- length(y)
weights <- tapply(weights, y, sum)
   cumu <- cumsum(weights)
   ecdf <- (cumu)/(cumu[length(cumu)]+1)
    fun <- stepfun(ysort, c(0, ecdf))
  class(fun) <- c("ecdf", "stepfun")
    return(fun)
  }
}
################################################################################
################################################################################
################################################################################
################################################################################
ESURV <- function (y, weights=NULL) 
{
  if (is.null(weights))
  {
    ysort <- unique(sort(y))
        n <- length(y)
     esurv <- 1-cumsum(tabulate(match(y, ysort)))/(n + 1)
       fun <- stepfun(ysort, c(1, esurv))
class(fun) <- c("esurv", "stepfun")
    return(fun)
  } else
  {
if (length(weights)!=length(y)) stop("y and weights have different length")
      i <- is.na(weights) | weights == 0
if (any(i)) 
    {
            y <- y[!i]
      weights <- weights[!i]
    } 
        ysort <- unique(sort(y))
            n <- length(y)       
      weights <- tapply(weights, y, sum)
         cumu <- cumsum(weights)
         esurv <- 1-(cumu)/(cumu[length(cumu)]+1)
    fun <- stepfun(ysort, c(1, esurv))
    class(fun) <- c("esurv", "stepfun")
    return(fun)
  }   
}  
# --------------------------------------------------------------
# this function is not exported 
#---------------------------------------------------------------
################################################################################
################################################################################
################################################################################
################################################################################
# This function allow weighted quantiles
# It is bsed on the function wtd.quantile() of package Hmisc
quantile_weights <- function (y, weights = NULL, probs = c(0, 0.25, 0.5, 0.75, 1))
{
  if (!length(weights))
    return(quantile(y, probs = probs))
  if (any(probs < 0 | probs > 1))
    stop("Probabilities must be between 0 and 1 inclusive")
  i <- is.na(weights) | weights == 0
  if (any(i))
  {
    y <- y[!i]
    weights <- weights[!i]
  }
  ysort <- unique(sort(y))
  weights1 <- tapply(weights, y, sum)
  cumu <- cumsum(weights1)
  x <- ysort
  wts <- weights1
  n <- sum(wts)
  order <- 1 + (n - 1) * probs
  low <- pmax(floor(order), 1)
  high <- pmin(low + 1, n)
  order <- order%%1
  allq <- approx(cumsum(wts), x, xout = c(low, high), method = "constant",
                 f = 1, rule = 2)$y
  k <- length(probs)
  quantiles <- (1 - order) * allq[1:k] + order * allq[-(1:k)]
  nams <- paste(format(round(probs * 100,
            if (length(probs) > 1) 2 - log10(diff(range(probs))) else 2)),
            "%", sep = "")
  names(quantiles) <- nams
  return(quantiles)
}
################################################################################
################################################################################
################################################################################
################################################################################
