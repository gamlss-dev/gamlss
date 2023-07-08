############################################################################
############################################################################
############################################################################
# the bucket plot
# TO DO
############################################################################
############################################################################
############################################################################  
# Fernanda De Bastiani, Mikis Stasinopoulos and Bob Rigby  
############################################################################
############################################################################
############################################################################ 
bp <- function (      obj = NULL, # this can be a variable or gamlss model or a list of both types  
                  weights = NULL, # this  needs to be checked     
                     type = c("moment", "centile.central", "centile.tail"), 
		                 xvar = NULL, # two version simple or formula
                bootstrap = TRUE, # whether to plot it
             no.bootstrap = 99,   # the number of bootstrap samples   
            col.bootstrap = c("lightblue", "pink", "khaki","thistle", "tan", "sienna1", "steelblue", "coral", "gold", "cyan"), # taken fro colors()  
            pch.bootstrap = rep(21,10), # design up to 10  the point type of bootstrap samples
              asCharacter = TRUE, # whether to plot name or point 
                col.point = rep("black", 10), # the colour of the plotted point or text
                pch.point = 1:10, #
                lwd.point = 2, #  
             text.to.show = NULL, # which text to show, the character of the y or something else 
                 cex.text = 1.5,  # size of text
                 col.text = "black", # color of text
              show.legend = FALSE,  # adding the legend this do not apply for multiple bp
		              n.inter = 4, # numner of cuts only for multiple bp
		          xcut.points = NULL, # if providing cuts only for multiple bp
		              overlap = 0, # if you overlaping the intervals only for multiple bp
		           show.given = TRUE, # only for multiple bp whether to show bars
		                  cex = 1, # pass it to  coplot()
		                  pch = 21, # pass it to  coplot()
                     data = NULL, 
                   bar.bg = c(num ="lightblue", fac="pink"),# colour or bars
		                  ...) # argument for coplot()
{
#############################################################################
#############################################################################
# local functions
#   i)   momBucket():  to draw the moment bucket used by panel.fun()
#  ii)   draw.Jarque.Bera(): used by panel.fun() 
# iii)   CI95() used by draw.Jarque.Bera()  
#  iv)   panel.fun():   used by bp() in coplot()
#   v)   check.overlap():  used by bp()
#  vi)   get.intervals(): used by bp()
# vii)   deparen(): used by bp()
##############################################################################
##############################################################################
# FUNCTION 1 
# This function  is used by panel.fun() to draw the moment bucket  
momBucket <- function(x)
{
         n <- length(x)
  boundary <- function(lty = 1, lwd = 2, col = 1){}
   fEGB2_1 <- function(x, deriv = 0L){}
   fEGB2_2 <- function(v){}
      fJSU <- function(x, deriv = 0L){}
     fSEP3 <- function(x, deriv = 0L){}
    fST3_1 <- function(x, deriv = 0L){}
    fST3_2 <- function(x, deriv = 0L){}
   fSHASHo <- function(x, deriv = 0L){}
load(system.file("doc", "MomentSkewKurt1.RData", package="gamlss.dist"))  
# logNO-JSU
  curve(fJSU, 0,1, ylim=c(-1,1), xlim=c(-1,1), xlab="transformed moment skewness",
        ylab="transformed moment kurtosis", col=gray(.7), lwd=2, lty=2, add=TRUE)
  flipfJSU <- function(x) fJSU(-x)
  curve( flipfJSU, -1,0, add=TRUE, col=gray(.7), lwd=2, lty=2)
# SHASHo
  curve(fSHASHo, 0.001,1, lty=3, col=gray(.4), lwd=2, add=TRUE)
  flipfSHASHo <- function(x) fSHASHo(-x)
  curve( flipfSHASHo, -1,-0.001, add=TRUE, lty=3, col=gray(.4), lwd=2)
# SEP3
  curve(fSEP3, 0.01,1, lty=4, col=gray(.6), lwd=2, add=TRUE)
  flipfSEP3 <- function(x) fSEP3(-x)
  curve(flipfSEP3, -1,-0.01, lty=4, col=gray(.6), lwd=2, add=TRUE)
# ST3
  curve(fST3_2, 0.5,1, lty=5, col=gray(.8), lwd=2,add=TRUE)
  flipfST3_2 <- function(x) fST3_2(-x)
  curve(flipfST3_2, -1,-0.5, lty=5, col=gray(.8), lwd=2, add=TRUE)
  curve(fST3_1, 0,.49, add=TRUE,  lty=5, col=gray(.8), lwd=2)
  flipfST3_1 <- function(x) fST3_1(-x)
  curve(flipfST3_1, -0.49,0,  lty=5, col=gray(.8), lwd=2, add=TRUE)
# EGB2
  curve(fEGB2_1, 0, .6666, lty=6, col=gray(.5), lwd=2, add=TRUE)
  flipfEGB2_1 <- function(x) fEGB2_1(-x)
  curve(flipfEGB2_1, -0.6666, 0, lty=6, col=gray(.5), lwd=2, add=TRUE)
  curve(fEGB2_2 , 0.0, .6666, lty=6, col=gray(.5), lwd=2, add=TRUE)
  flipfEGB2_2 <- function(x) fEGB2_2(-x)
  curve(flipfEGB2_2, -0.6666, 0,  lty=6, col=gray(.5), lwd=2, add=TRUE)
# local
  boundary()
# local function
boundaryR<-function()
  {
    tskew<-seq(0,0.99999,length=101)
    # tskew<-seq(-0.99999,0,length=1000)
     skew <- tskew/(1-tskew)
     kurt <- 1 + (skew^2)
    ekurt <- kurt - 3
    tkurt <- ekurt/(1+abs(ekurt))
lines(tkurt~I(-tskew), type="l", lty=1, lwd=2, col=1)
  }
  boundaryR()
  grid()
  lines(c(-1,1), c(1,1), lwd=2)
  lines(c(-1,1), c(1,1), lwd=2)  
}
##############################################################################
##############################################################################
# FUNCTION 2
# This function  is used by panel.fun() to draw the boundaries of Jarque.Bera 
draw.Jarque.Bera <- function(n) 
{
  xx <- seq(0, 0.99, length = 101)
  y1 <- CI95(xx, n=n)
  y2 <- CI95.2(-xx, n=n)
  y3 <- CI95.3(-xx, n=n)
  y4 <- CI95.4(xx, n=n)
  xy <- na.omit(data.frame(x = c(xx, xx[101:1], -xx, -xx[101:1]), 
                           y = c(y1, y4[101:1], y3, y2[101:1])))
  polygon(xy$x, xy$y, lty = 11, col = gray(0.97))  
}
##############################################################################
##############################################################################
# FUNCTION 3
# used by draw.Jarque.Bera
CI95 <- function(x,n) {
  if (any(abs(x) >= 1)) 
    stop(" x should be in (-1,1)")
  gamma.1 <- x/(1 - x)
  gamma.2 <- ((24/n) * (qchisq(0.95, df = 2) - (n/6) * 
                          gamma.1^2))^0.5
  gamma.2t <- gamma.2/(1 + abs(gamma.2))
  gamma.2t
}
CI95.2 <- function(x,n) CI95(-x,n)
CI95.3 <- function(x,n) -CI95(-x,n)
CI95.4 <- function(x,n) -CI95(x,n)  
##############################################################################
##############################################################################
# FUNCTION 4 This is the important function
panel.fun <- function (x, y, ...) 
	{
   id <- get("id", envir = parent.frame(1))
    n <- length(y)  
#    N <- length(id)
#    if ((n/N)*100 <1) warning("only ", n, " out of ", N, " observations are in this panel \n")
# if moment ##################################################################    
if (type=="moment")
{
  momBucket(Y[[1]]) # draw the bucket
  draw.Jarque.Bera(n)
for (i in 1:K)
 {
  y  <- Y[[i]][id]
  w  <- W[[i]][id]  
if (bootstrap)
  {
    for (j in 1:no.bootstrap)
    {
      ind <- sample.int(n, n,  replace = TRUE)
      sk <-  momentSK(y[ind], weights=w[ind])
      points(sk$trans.mom.skew, sk$trans.mom.kurt, pch=pch.bootstrap[i] , 
             col=col.bootstrap[i])
    }
  }
}
# this loop is repeated here for the symbols to be clear   
  for (i in 1:K)
  {
    y  <- Y[[i]][id]
    w  <- W[[i]][id]  
    sk <-  momentSK(y, weights=w)
  if  (asCharacter)
    text(sk$trans.mom.skew, sk$trans.mom.kurt, text.to.show[i], cex=cex.text, col=col.text ) 
    else 
    points(sk$trans.mom.skew, sk$trans.mom.kurt, col=col.point[i], pch=pch.point[i],
         lwd=lwd.point) 
    text(-0.9,-0.9, n)
  }  
}
# end if moment ##############################################################
# if centile central##########################################################       
if (type=="centile.central")
{
  id <- get("id", envir = parent.frame(1))
   n <- length(y)  
  cenCentralBucket() # draw the bucket
for (i in 1:K)
  {  
  y  <- Y[[i]][id]
  w  <- W[[i]][id]
    if (bootstrap)
   {
     for (j in 1:no.bootstrap)
     {
       ind <- sample.int(n, n,  replace = TRUE)
       sk <-  centileSK(y[ind], weights=w[ind])
      points(sk$S0.25,  sk$trans.K0.01,  pch=pch.bootstrap[i], col=col.bootstrap[i])
     }
    }
} 
# this loop is repeated here for the symbols to be clear   
for (i in 1:K)
  {
    y  <- Y[[i]][id]
    w  <- W[[i]][id]  
    sk <-centileSK(y, weights=w)
    text(-0.9,-0.9, n)
     if  (asCharacter)
        text(sk$S0.25, sk$trans.K0.01, text.to.show[i], cex=cex.text, col=col.text ) 
      else 
        points(sk$S0.25, sk$trans.K0.01, col=col.point[i], pch=pch.point[i],
               lwd=lwd.point) 
  } 
  
  }    
# end centile central#########################################################    
# if centile tailt ###########################################################       
if (type=="centile.tail")
{
  id <- get("id", envir = parent.frame(1))
   n <- length(y)  
  cenTailBucket()  # draw the bucket
for (i in 1:K)
  {  
   y  <- Y[[i]][id]
   w  <- W[[i]][id]
    if (bootstrap)
    {
      for (j in 1:no.bootstrap)
      {
        ind <- sample.int(n, n,  replace = TRUE)
        sk <-  centileSK(y[ind], weights=w[ind])
        points(sk$S0.01, sk$trans.K0.01,  pch=pch.bootstrap[i], col=col.bootstrap[i])
      }
    }
} 
  for (i in 1:K)
  {
    y  <- Y[[i]][id]
    w  <- W[[i]][id]  
    sk <-centileSK(y, weights=w)
    text(-0.9,-0.9, n)
    if  (asCharacter)
      text(sk$S0.01, sk$trans.K0.01, text.to.show[i], cex=cex.text, col=col.text ) 
    else 
      points(sk$S0.01, sk$trans.K0.01, col=col.point[i], pch=pch.point[i],
             lwd=lwd.point) 
  }    
}
    
}    
# if centile tail ###############################################  
##############################################################################
##############################################################################
# FUNCTION 5
check.overlap <- function(interval)
	{
		if (!is.matrix(interval) ) {stop(paste("The interval specified is not a matrix."))}
		if (dim(interval)[2] !=2) {stop(paste("The interval specified is not a valid matrix.\nThe number of columns should be equal to 2."))}
		crows = dim(interval)[1]
		for (i in 1:(crows-1))
		{
			#if (interval[i,2] != interval[i+1,1]) {interval[i+1,1]=interval[i,2]}
			if (!(abs(interval[i,2]-interval[i+1,1])<0.0001)) {interval[i+1,1]=interval[i,2]}
		}
		return(interval)
	}
##############################################################################
##############################################################################
# FUNCTION 6
# the problem here is that coplot uses for given.values      
#  MS Tuesday, September 21, 2004 
#  min <=  x <= 20                   min <=  x < 20             min <=  x <= 20-extra  
#   20 <=  x <= 30    instead   of    20 <=  x < 30  here uses   20 <=  x <= 30-extra  
#   30 <=  x <= max                   30 <=  x <= max            30 <=  x <= max+extra 
get.intervals <- function (xvar, xcut.points ) 
	{
		if (!is.vector(xcut.points))  {stop(paste("The interval is not a vector."))}
		if ( any((xcut.points < min(xvar)) | any(xcut.points > max(xvar))))
		{stop(paste("The specified `xcut.points' are not within the range of the x: (", min(xvar),
							" , ", max(xvar), ")"))}
		extra <- (max(xvar)-min(xvar))/100000
		  int <- c(min(xvar), xcut.points, (max(xvar)+2*extra))
		   ii <- 1:(length(int)-1)
		    r <- 2:length(int)
		   x1 <- int[ii]
		   xr <- int[r]-extra
		if (any(x1>xr)) {stop(paste("The interval is are not in a increasing order."))}
		cbind(x1,xr)
	}
################################################################################
################################################################################
# FUNCTION 7
# this function comes from coplot() to help to get the variables used if more that one xvar is used 
deparen <- function(expr) 
{
  while (is.language(expr) && !is.name(expr) && deparse(expr[[1L]])[1L] == 
           "(") expr <- expr[[2L]]
  expr
}
###############################################################################
###############################################################################
# FUNCTION 8
cenCentralBucket <- function()
{
  cEGB2_1 <- cEGB2_2 <-tEGB2_1 <-tEGB2_2 <- cJSU <-  tJSU <- 
   tST3_1 <- tST3_2 <-function (x, deriv = 0L) {}
  cSB <- tSB <- cSEP3  <- tSEP3 <- cSHASH <- tSHASH <- cST3_1 <- 
  cST3_2 <- function (x, deriv = 0L) {}
  cEGB2_1_data <- cEGB2_2_data <- data.frame( cskew=0, ckurt=0)
  load(system.file("doc", "CentileSkewKurt.RData", package="gamlss.dist"))     
    curve(cSEP3, 0.01,.99,  ylim=c(-1,1),  xlim=c(-1,1), xlab="central centile skewness",
          ylab="transformed centle kurtosis",
          lty=4, col=gray(.6), lwd=2, add=TRUE)
    flipcSEP3 <- function(x) cSEP3(-x)
    curve( flipcSEP3, -1,-0.01, add=TRUE, lty=4, col=gray(.6), lwd=2)
    curve(cST3_2, 0.1445, 1, add=T, lty=5, col=gray(.8), lwd=2)
    curve(cST3_1, 0, 0.1441, add=T, lty=5, col=gray(.8), lwd=2)
    flipcST3_2 <- function(x) cST3_2(-x)
    curve(flipcST3_2, -1, -0.144, add=T, lty=5, col=gray(.8), lwd=2)
    flipcST3_1 <- function(x) cST3_1(-x)
    curve(flipcST3_1, -0.1441,0,  add=T, lty=5, col=gray(.8), lwd=2)
    curve(cJSU, 0.01, .99, add=T, lty=2, col=gray(.7),  lwd=2)
    flipcJSU <- function(x) cJSU(-x)
    curve(flipcJSU, -.99, -0.01, add=T, lty=2, col=gray(.7),  lwd=2)
    curve(cSHASH, 0, .965, add=T, lty=3, col=gray(.4), lwd=2)
    flipcSHASH <- function(x) cSHASH(-x)
    curve(flipcSHASH, -.965, 0, add=T, lty=3, col=gray(.4), lwd=2)
    #lines(cEGB2_1Data$cskew,cEGB2_1Data$ckurt,  lty=6, col=gray(.5), lwd=2)
    #lines(-cEGB2_1Data$cskew,cEGB2_1Data$ckurt,  lty=6, col=gray(.5), lwd=2)
    lines(cEGB2_1_data$cskew,cEGB2_1_data$ckurt,  lty=6, col=gray(.5), lwd=2)
    lines(-cEGB2_1_data$cskew,cEGB2_1_data$ckurt,  lty=6, col=gray(.5), lwd=2)
    lines(cEGB2_2_data$cskew[c(-1,-2)],cEGB2_2_data$ckurt[c(-1,-2)],  lty=6, col=gray(.5), lwd=2) 
    lines(-cEGB2_2_data$cskew[c(-1,-2)],cEGB2_2_data$ckurt[c(-1,-2)],  lty=6, col=gray(.5), lwd=2)
    # curve(cEGB2_2, 0,  0.26, add=T,lty=6, col=gray(.5), lwd=2)
    # flipcEGB2_2 <- function(x) cEGB2_2(-x)
    # curve(flipcEGB2_2, -0.26, 0, add=T,lty=6, col=gray(.5), lwd=2)
    # lines(c(-0.26,-0.26), c(0.5078402, 0.5754391), col=gray(.5), lty=6,  lwd=2)
    curve(cSB, 0, 0.99, add=T, lty=7, col=gray(.9), lwd=2)
    flipcSB <- function(x) cSB(-x)
    curve(flipcSB, -0.99, 0, add=T, lty=7, col=gray(.9), lwd=2)
    lines(c(-1,1), c(1,1), lwd=2)
    grid()  
}
###############################################################################
###############################################################################
# FUNCTION 9
cenTailBucket <- function()
{
  cEGB2_1 <- cEGB2_2 <-tEGB2_1 <-tEGB2_2 <- cJSU <-  tJSU <- 
    tST3_1 <- tST3_2 <-function (x, deriv = 0L) {}
  cSB <- tSB <- cSEP3  <- tSEP3 <- cSHASH <- tSHASH <- cST3_1 <- 
    cST3_2 <- function (x, deriv = 0L) {}
  cEGB2_1_data <- cEGB2_2_data <- data.frame( cskew=0, ckurt=0)
  load(system.file("doc", "CentileSkewKurt.RData", package="gamlss.dist"))   
  curve(tSEP3, 0.01,.99,  ylim=c(-1,1),  xlim=c(-1,1), xlab="tail centile skewness",
        ylab="transformed centle kurtosis",
        lty=4, col=gray(.6), lwd=2, add=TRUE)
  fliptSEP3 <- function(x) tSEP3(-x)
  curve( fliptSEP3, -1,-0.01, add=TRUE, lty=4, col=gray(.6), lwd=2)
  curve(tSHASH, 0, .996, add=TRUE, lty=3, col=gray(.4), lwd=2)
  fliptSHASH <- function(x) tSHASH(-x)
  curve(fliptSHASH, -.996,0, add=T, lty=3, col=gray(.4), lwd=2)
  curve(tST3_2, 0.484, 1, add=TRUE, lty=5, col=gray(.8), lwd=2)
  fliptST3_2 <- function(x) tST3_2(-x)
  curve(fliptST3_2, -1,  -0.484,  add=T, lty=5, col=gray(.8), lwd=2)
  curve(tST3_1,  0, 0.482, add=TRUE, lty=5, col=gray(.8), lwd=2)
  fliptST3_1 <- function(x) tST3_1(-x)
  curve(fliptST3_1, -0.482, 0,  add=T, lty=5, col=gray(.8), lwd=2)
  curve(tJSU, 0.01 ,0.99, add=T, lty=2, col=gray(.7),  lwd=2)
  fliptJSU <- function(x) tJSU(-x)
  curve(fliptJSU,  -0.99, -0.01, add=T, lty=2, col=gray(.7),  lwd=2)
  curve(tEGB2_1, 0,  0.70, add=T, lty=6, col=gray(.5), lwd=2)
  fliptEGB2_1 <- function(x) tEGB2_1(-x)
  curve(fliptEGB2_1, -0.70, 0, add=T, lty=6, col=gray(.5), lwd=2)
  curve(tEGB2_2, 0,  0.70, add=T, lty=6, col=gray(.5), lwd=2)
  fliptEGB2_2 <- function(x) tEGB2_2(-x)
  curve(fliptEGB2_2, -0.70, 0, add=T, lty=6, col=gray(.5), lwd=2)
  curve(tSB, 0,1, add=T, lty=7, col=gray(.9), lwd=2)
  fliptSB <- function(x) tSB(-x)
  curve(fliptSB, -1,0, add=T, lty=7, col=gray(.9), lwd=2)
  grid()
  lines(c(-1,1), c(1,1), lwd=2) 
}
# enf of lcal functions  
###############################################################################
###############################################################################
###############################################################################
###############################################################################
##-----------------------------------------------------------------------------
##-here is the main function --------------------------------------------------
##-----------------------------------------------------------------------------
###############################################################################
###############################################################################
###############################################################################
###############################################################################
type <-  match.arg(type)
## get residuals
# check whether a list
#DataExist <- FALSE
        K <- 1L
        Y <- n <- W <-list()
if (is(obj,"list"))
{
  K <- length(obj)
  for (i in 1L:K)
  {
    Y[[i]] <-   if (is(obj[[i]],"gamlss")) resid(obj[[i]])
    else obj[[i]]
    W[[i]] <-   if (is(obj[[i]],"gamlss")) obj[[i]]$weight
    else if (is.null(weights)) rep(1, length(Y[[i]])) else weights[[i]]
    n[[i]] <-   length(Y[[i]])
     NAMES <-  paste(substitute(obj))[-1]
  }
} else
{
  Y[[1L]] <-  if (is(obj,"gamlss")) resid(obj)
  else obj  
  W[[1L]] <-  if (is(obj,"gamlss")) obj$weights
  else if (is.null(weights)) rep(1L, length(obj)) else weights
  n[1L] <- length(Y[[1L]])
  NAMES <-  paste(substitute(obj))
}  
text.to.show <- if (is.null(text.to.show)) NAMES else  text.to.show
##----------------------------------------------------------------------------
##############################################################################
##############################################################################
##----------------------------------------------------------------------------
# case 1 no x-variable then just somple plot
if(is.null(xvar)) # if xvar=NULL
	{
  for (i in 1L:K)
  {
if (i==1L) add <- FALSE else  add <- TRUE
#if (asCharacter==TRUE) 
#  text.to.show <-if (K>1L) paste(substitute(obj)[[1L+i]]) else paste(substitute(obj))  
if (type=="moment")
    {
    out <-checkMomentSK(Y[[i]], weights=W[[i]],
              add = add, bootstrap = bootstrap, no.bootstrap = no.bootstrap, 
              col.bootstrap = col.bootstrap[i], pch.bootstrap = pch.bootstrap[i], 
              asCharacter = asCharacter, col.point = col.point[i] , pch.point = pch.point[i], 
              lwd.point = lwd.point,  text.to.show = text.to.show[i], cex.text = cex.text, 
              col.text = col.text, show.legend = show.legend)
    }
if (type=="centile.central")
    {
      out <-checkCentileSK(Y[[i]], weights=W[[i]], type='central',
              add = add, bootstrap = bootstrap, no.bootstrap = no.bootstrap, 
              col.bootstrap = col.bootstrap[i], pch.bootstrap = pch.bootstrap[i], 
              asCharacter = asCharacter, col.point = col.point[i] , pch.point = pch.point[i], 
              lwd.point = lwd.point,   text.to.show = text.to.show[i], cex.text = cex.text, 
              col.text = col.text, show.legend = show.legend)
    }
if (type=="centile.tail")
    {
      out <-checkCentileSK(Y[[i]], weights=W[[i]], type='tail',
              add = add, bootstrap = bootstrap, no.bootstrap = no.bootstrap, 
              col.bootstrap = col.bootstrap[i], pch.bootstrap = pch.bootstrap[i], 
              asCharacter = asCharacter, col.point = col.point[i] , pch.point = pch.point[i], 
              lwd.point = lwd.point,   text.to.show = text.to.show[i], cex.text = cex.text, 
              col.text = col.text, show.legend = show.legend)
    } 
  }
  return(invisible(out))   
}  
##----------------------------------------------------------------------------
##############################################################################
##############################################################################
##----------------------------------------------------------------------------
## get data if exist
if (!is.null(data)) Data <- data
  form <- if (is(xvar,"formula"))   {as.formula(paste0("y.y~x.x|", as.character(xvar))[2])}
         else {as.formula(paste0("y.y~x.x|", deparse(substitute(xvar))))}
## get weights
   y.y <- Y[[1L]]   
   x.x <- W[[1L]]  
## if w=0 reduce if w=freq expand if all(w=1L) same
## here we need the xvar   
#if (all(trunc(w)==w))  xvar <- rep(xvar, w) 
if (!is(xvar,"formula"))
{
  if(is.null(xcut.points)) 
  { # getting the intervals automatic
    if (is(xvar,"factor")) stop("factors should be used with formula, i.e. ~f1")
    given.in <- co.intervals(xvar, number=n.inter, overlap=overlap )
    if (overlap==0) given.in <- check.overlap(given.in) 
  }                 
  else
  { # if xcut.points is set
    given.in <- get.intervals(xvar, xcut.points )
  }
  coplot(form, 
         given.values = given.in, 
         panel = panel.fun,
         ylim = c(-1, 1),
         xlim = c(-1, 1),
         ylab = "skewness",
         xlab = "kurtosis",
         show.given = show.given,
         pch = pch,
         cex = cex, 
         bar.bg =  bar.bg ,...)
  return(invisible())
  } else
  {
    coplot(form, 
   # given.values = given.in, 
           panel = panel.fun,
           ylim = c(-1, 1),
           xlim = c(-1, 1),
           ylab = "skewness",
           xlab = "kurtosis",
           show.given = show.given,
           pch = pch,
           cex = cex, 
           number = n.inter,
           bar.bg =  bar.bg ,...)
    return(invisible())
  }
 
}#--------------------------------------------------------- END FUNCTION
#############################################################################
#############################################################################
#############################################################################