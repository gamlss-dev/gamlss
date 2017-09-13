## this is a  cycle Penalized Beta  regression splines smoother
## Paul Eilers, Mikis Stasinopoulos and Bob Rigby
## last modified Saturday, August 28, 2009 
#----------------------------------------------------------------------------------------
pbc<-function(x, df = NULL, lambda = NULL, control=pbc.control(...), ...) 
{
## this function is based on Paul Eilers' penalised beta regression splines function
## lambda : is the smoothing parameter
## df : are the effective df's
## if both lambda=NULL  and df=NULL then lambda is estimated using the different method
## methods are "ML", "ML-1", "EM", "GAIC" and "GCV"  
## if df is set to number but lambda is NULL then df are used for smoothing
## if lambda is set to a number (whether df=NULL  or not) lambda is used for smoothing
# ---------------------------------------------------
## local function
## creates the basis for p-splines
## Paul Eilers' function
#-------------------------------------------------------------------------------
bbase <- function(x, xl = min(x), xr = max(x), nseg = 10, deg = 3)
{
  tpower <- function(x, t, p)
    # Truncated p-th power function
    (x - t) ^ p * (x > t)
  # Construct B-spline basis
  dx <- (xr - xl) / nseg
  knots <- seq(xl - deg * dx, xr + deg * dx, by = dx)
  P <- outer(x, knots, tpower, deg)
  n <- dim(P)[2]
  D <- diff(diag(n), diff = deg + 1) / (gamma(deg + 1) * dx ^ deg)
  B <- (-1) ^ (deg + 1) * P %*% t(D)
  B
}
#-------------------------------------------------------------------------------
## Paul Eilers' function
cbase <- function(x, xl, xr, ndx, deg)
{
  # Construct circular B-spline basis
  # Domain: xl to xr, number of segmants on domain: ndx,  degree: deg
  # Wrap around to cyclic basis
  B0 = bbase(x, xl = xl, xr = xr, nseg = ndx, deg = deg)
  n = ncol(B0) - deg # 13-3=10
  cc = (1:deg) + n   # 11 12 13
  B = B0[, 1:n]      # B is 100 X 10 
  B[, 1:deg] = B[, 1:deg] + B0[, cc] 
  # B[,c(1,2,3)] = B[,c(1,2,3)] + B[,c(13,12,11)] Friday, BR and MS October 9, 2009 
  B
}
#------------------------------------------------------------------------------
## Paul Eilers' function
## Modified by Mikis
## the Cdiff function allows for order= 1 and 2 2
## for order 2 allows ti have sinoide behaviour 
## 2* cos(2 * pi / n or not
Cdiff <- function(n, order=2, sin=TRUE)
{
  if (order==1)
  {
    D2 <- matrix(0, n, n + 1)
    p <- c(1, -1)
    for (k in 1:n) D2[k, c(0:1) + k] = p
    D <- D2[, 2:(n + 1)]
    D[, n] <- D[, n] + D2[, 1]
  }
  if (order==2)
  {
    D2 <- matrix(0, n, n + 2)
    p <- if (sin==TRUE)  c(-1, 2* cos(2 * pi / n), -1) else c(-1, 2, -1)
    for (k in 1:n) D2[k, (0:2) + k] = p
    D <- D2[, 2:(n + 1)]
    D[, 1] <- D[, 1] + D2[, n + 2]
    D[, n] <- D[, n] + D2[, 1]
  }
  D
}
#-------------------------------------------------------------------------------
#--------------------------------------------------------------------------------
# the main function starts here
         scall <- deparse(sys.call())
            lx <- length(x)
if (is(x,"factor")) # ||control$ts==TRUE
 {
          xval <- as.numeric(unique(x))
            xl <- min(xval)
            xr <- max(xval)
            nl <- nlevels(x) 
          # function(x, xl, xr, ndx, deg)
             X <- cbase(as.numeric(x), xl, xr,  nl, 1)  # create the basis
 }            
else
 {
   no.dist.val <-  length(table(x))
 control$inter <- if (lx<100) 10 else control$inter 
 control$inter <- if (no.dist.val<=control$inter) no.dist.val else control$inter 
            xl <- min(x)
            xr <- max(x)
          xmax <- xr #+ 0.01 * (xr - xl) # BR and MS Friday, October 9, 2009 
          xmin <- xl #- 0.01 * (xr - xl)
             X <- cbase(x, xl, xr, control$inter, control$degree)  # create the basis
 }
          # Cyclic penalty
            nb <- ncol(X)
             D <- Cdiff(nb, order=control$order, sin=control$sin)
            # D <- diff(D)  # not know yt if
             if(!is.null(df)) # degrees of freedom
             {
             if (df>(dim(X)[2]-2)) 
              {df <- 3;  
              warning("The df's exceed the number of columns of the design matrix", "\n",  "   they are set to 3") }
              df <- if (df < 1)  1  else  df+1
              if (df < 1)  warning("the df are set to 1")    
             }
## here we get the gamlss environment and a random name to save
## the starting values for lambda within gamlss()
## get gamlss environment
#--------
     rexpr<-regexpr("gamlss",sys.calls())
for (i in 1:length(rexpr)){ 
    position <- i 
    if (rexpr[i]==1) break}
gamlss.environment <- sys.frame(position)
#--------
## get a random name to use it in the gamlss() environment
#--------
               sl <-sample(letters, 4)
      fourLetters <- paste(paste(paste(sl[1], sl[2], sep=""), sl[3], sep=""),sl[4], sep="")
  startLambdaName <- paste("start.Lambda",fourLetters, sep=".")
## put the starting values in the gamlss()environment
#--------
   assign(startLambdaName, control$start, envir=gamlss.environment)
#--------
          xvar <- rep(0,length(x)) # 
      attr(xvar, "control")       <- control
      attr(xvar, "D")             <- D
      attr(xvar, "X")             <- X
      attr(xvar, "df")            <- df 
      attr(xvar, "call")          <- substitute(gamlss.pbc(data[[scall]], z, w)) 
      attr(xvar, "lambda")        <- lambda
      attr(xvar, "gamlss.env")    <- gamlss.environment
      attr(xvar, "NameForLambda") <- startLambdaName
      attr(xvar, "class")         <- "smooth"
      xvar
}
#----------------------------------------------------------------------------------------
# control function for cy()
##---------------------------------------------------------------------------------------
pbc.control <- function(inter = 20, degree= 3, order = 2, start=10, 
                       method=c("ML","GAIC", "GCV"), k=2, sin=TRUE, ...)
{ 
##  Control function for cy()
##  MS  Tuesday, March 24, 2009
## inter : is the number of equal space intervals in x  
## degree: is the degree of the polynomial 
## order refers to differences in the penalty for the coeficients 
## order = 0 : white noise random effects
## order = 1 : random walk
## order = 2 : random walk of order 2
## order = 3 : random walk of order 3
        if(inter <= 0) {
warning("the value of inter supplied is less than 0, the value of 10 was used instead")
                inter <- 10 }
        if(degree <= 0) {
warning("the value of degree supplied is less than zero or negative the default value of 3 was used instead")
                degree <- 3}                
        if(order < 0) {
warning("the value of order supplied is zero or negative the default value of 2 was used instead")
                order <- 2}
       if(order >= 3) {
  warning("the value of order supplied is greater than 2 the  default value of 2 was used instead")
  order <- 2}
        if(k <= 0) {
warning("the value of GAIC/GCV penalty supplied is less than zero the default value of 2 was used instead")
                k <- 2}   
method <- match.arg(method)                          
        list(inter = inter, degree = degree,  order = order, start=start, method= method, k=k,  sin=as.logical(sin)[1])
}
#----------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------
gamlss.pbc <- function(x, y, w, xeval = NULL, ...)
{
# -------------------------------------------------- 
# functions within
# a siple penalized regression
regpen <- function(y, X, w, lambda, D)
  {
             G <- lambda * t(D) %*% D
            XW <- w * X
           XWX <- t(XW) %*% X
          beta <- solve(XWX + G, t(XW) %*% y)
            fv <- X %*%beta
             H <- solve(XWX + G, XWX)
         #  edf <- sum(diag(H))
           fit <- list(beta = beta, edf = sum(diag(H)))
  return(fit)  
  }
#-------------------------------------------------------------------------------
regpen <- function(y, X, w, lambda, D)# original
{
  #         p <- dim(X)[2]
  #       qrX <- qr(sqrt(w)*X, tol=.Machine$double.eps^.8)  
  #         R <- qr.R(qrX)
  RD <- rbind(R,sqrt(lambda)*D) # 2p x p matrix 
  svdRD <- svd(RD)                 # U 2pxp D pxp V pxp
  ##             take only the important values    
  rank <- sum(svdRD$d>max(svdRD$d)*.Machine$double.eps^.8)
  U1 <- svdRD$u[1:p,1:rank]     # U1 p x rank 
  # I am not sure what are consequances in introducing this ???
  y1 <- t(U1)%*%Qy #  t(Q)%*%(sqrt(w)*y)        # rankxp pxn nx1 => rank x 1 vector 
  #     beta <- svdRD$v[,1:rank] %*%diag(1/svdRD$d[1:rank])%*%y1 
  beta <- svdRD$v[,1:rank] %*%(y1/svdRD$d[1:rank])
  #        1/(svdRD$d^2)
  #print((svdRD$v)%*%t(svdRD$v), digits=1)
  HH <- (svdRD$u)[1:p,1:rank]%*%t(svdRD$u[1:p,1:rank])
  df <- sum(diag(HH))
  fit <- list(beta = beta, edf = df)
  return(fit)  
}
#-------------------------------------------------------------------------------
#--------------------------------------------------
## function to find lambdas miimizing the local GAIC        
     fnGAIC <- function(lambda, k)
    {
       fit <- regpen(y=y, X=X, w=w, lambda=lambda, D)
        fv <- X %*% fit$beta         
      GAIC <- sum(w*(y-fv)^2)+k*fit$edf 
    # cat("GAIC", GAIC, "\n")
      GAIC   
    }
#-------------------------------------------------------------------------------
# ## function to find lambdas miimizing the local GAIC        
fnGAIC <- function(lambda, k)
{
  fit <- regpen(y=y, X=X, w=w, lambda=lambda, D)
  fv <- X %*% fit$beta         
  GAIC <- sum(w*(y-fv)^2)+k*fit$edf 
  # cat("GAIC", GAIC, "\n")
  GAIC   
}
# 
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# ## function to find the lambdas which minimise the local GCV 
fnGCV <- function(lambda, k)
{
  I.lambda.D <- (1+lambda*UDU$values)
  edf <- sum(1/I.lambda.D)
  y_Hy2 <- y.y-2*sum((yy^2)/I.lambda.D)+sum((yy^2)/((I.lambda.D)^2))
  GCV <- (n*y_Hy2)/(n-k*edf)^2
  GCV
}  
#------------------------------------------------------------------------------ 
## local function to get df using eigen values
    edf1_df <- function(lambda)
           {
           edf <-  sum(1/(1+lambda*UDU$values))
           (edf-df)
           }  
#------------------------------------------------------------------
# the main function starts here
# get the attributes
         X <-  if (is.null(xeval)) as.matrix(attr(x,"X")) #the trick is for prediction
              else  as.matrix(attr(x,"X"))[seq(1,length(y)),]
         D <- as.matrix(attr(x,"D")) # penalty
    lambda <- as.vector(attr(x,"lambda")) # lambda
        df <- as.vector(attr(x,"df")) # degrees of freedom
   control <- as.list(attr(x, "control")) 
gamlss.env <- as.environment(attr(x, "gamlss.env"))
startLambdaName <- as.character(attr(x, "NameForLambda")) 
     order <- control$order # the order of the penalty matrix
         N <- sum(w!=0) # DS+FDB 3-2-14    
         n <- nrow(X) # the no of observations
         p <- ncol(D) # the rows of the penalty matrix
       qrX <- qr(sqrt(w)*X, tol=.Machine$double.eps^.8)  
        R <- qr.R(qrX)
        Q <- qr.Q(qrX) 
      Qy  <- t(Q)%*%(sqrt(w)*y)
     tau2  <- sig2 <- NULL
# now the action depends on the values of lambda and df
#-------------------------------------------------------------------- 
  lambdaS <-  get(startLambdaName, envir=gamlss.env) ## geting the starting value
 if (lambdaS>=1e+07) lambda <- 1e+07 # MS 19-4-12
 if (lambdaS<=1e-07) lambda <- 1e-07 # MS 19-4-12 
# case 1: if lambda is known just fit
 if (is.null(df)&&!is.null(lambda)||!is.null(df)&&!is.null(lambda))
 {
          fit <- regpen(y, X, w, lambda,  D)
           fv <- X %*% fit$beta        
 } # case 2: if lambda is estimated ------------------------------------------- 
 else if (is.null(df)&&is.null(lambda)) 
 { #   
  # cat("----------------------------","\n")
  lambda <- lambdaS  # MS 19-4-12## geting the starting value
  # if ML ----------------------
switch(control$method,
    "ML"={
      for (it in 1:50) 
           {
              fit  <- regpen(y, X, w, lambda, D) # fit model
            gamma. <- D %*% as.vector(fit$beta)  # get the gamma differences
                fv <- X %*% fit$beta             # fitted values
              sig2 <- sum(w * (y - fv) ^ 2) / (N - fit$edf) # DS+FDB 3-2-14
             tau2 <- sum(gamma. ^ 2) / (fit$edf-order)# see LNP page 279
             if(tau2<1e-7) tau2 <- 1.0e-7 # MS 19-4-12
        lambda.old <- lambda
            lambda <- sig2 / tau2 # maybe only 1/tau2 will do since it gives exactly the EM results see LM-1
            if (lambda<1.0e-7) lambda<-1.0e-7 # DS Saturday, April 11, 2009 at 14:18
            if (lambda>1.0e+7) lambda<-1.0e+7 # DS 29 3 2012
             #    cat("iter tau2 sig2",it,tau2, sig2, '\n')
            if (abs(lambda-lambda.old) < 1.0e-7||lambda>1.0e10) break
            assign(startLambdaName, lambda, envir=gamlss.env)
             #cat("lambda",lambda, '\n')
           }
         },
         "GAIC"=  #--------------------------------------------------------------- GAIC
   {
     lambda <- nlminb(lambda, fnGAIC,  lower = 1.0e-7, upper = 1.0e7, k=control$k)$par 
        fit <- regpen(y=y, X=X, w=w, lambda=lambda, D)
         fv <- X %*% fit$beta     
    assign(startLambdaName, lambda, envir=gamlss.env)
   },
"GCV"={   #-------------------------------------------------------------- GCV
  # 
     wy <- sqrt(w)*y
    y.y <- sum(wy^2)
   Rinv <- solve(R)
      S <- t(D)%*%D
    UDU <- eigen(t(Rinv)%*%S%*%Rinv)
     yy <- t(UDU$vectors)%*%Qy #t(qr.Q(QR))%*%wy
 lambda <- nlminb(lambda, fnGCV,  lower = 1.0e-7, upper = 1.0e7, k=control$k)$par
    fit <- regpen(y=y, X=X, w=w, lambda=lambda, D)
     fv <- X %*% fit$beta     
  assign(startLambdaName, lambda, envir=gamlss.env) 
})
 }
  else # case 3 : if df are required---------------------------------
  { 
    Rinv <- solve(R)
     S   <- t(D)%*%D
     UDU <- eigen(t(Rinv)%*%S%*%Rinv)           
  lambda <- if (sign(edf1_df(0))==sign(edf1_df(100000))) 100000  # in case they have the some sign
    else  uniroot(edf1_df, c(0,100000))$root
    # if (any(class(lambda)%in%"try-error")) {lambda<-100000}   
    fit <- regpen(y, X, w, lambda, D)
    fv <- X %*% fit$beta
  }#--------------------------------------------------------------------------end of case 3
  
  #  but this is working
#Version 4 -------------------------------------------------- 
   waug <- as.vector(c(w, rep(1,nrow(D))))
   xaug <- as.matrix(rbind(X,sqrt(lambda)*D))
    lev <- hat(sqrt(waug)*xaug,intercept=FALSE)[1:n] # get the hat matrix
# the variance of the smootherz
     # se <-  sqrt(diag(solve(XWX + lambda * t(D) %*% D)))
coefSmo <- list(   coef = fit$beta,
                     fv = fv, 
                 lambda = lambda, 
                    edf = fit$edf, 
                  sigb2 = tau2, 
                  sige2 = sig2,
                   sigb = if (is.null(tau2)) NA else sqrt(tau2),
                   sige = if (is.null(sig2)) NA else sqrt(sig2),
                 method = control$method)
class(coefSmo) <- c("pbc", "pb") 
if (is.null(xeval)) # if no prediction 
    {
     list(fitted.values=fv, residuals=y-fv, var=lev, nl.df =fit$edf-1,
          lambda=lambda, coefSmo=coefSmo)
    }                            
else # for prediction 
    {
     ll <- dim(as.matrix(attr(x,"X")))[1]
     nx <- as.matrix(attr(x,"X"))[seq(length(y)+1,ll),]
   pred <- drop(nx %*% fit$beta)
   pred
    }    
}
#----------------------------------------------------------------------------------------
print.pbc  <- function (x, digits = max(3, getOption("digits") - 3), ...) 
{   
  cat("Cycle P-spline fit using the gamlss function cy() \n")
  cat("Degrees of Freedom for the fit :", x$edf, "\n")
  cat("Random effect parameter sigma_b:", format(signif(x$sigb)), "\n")  
  cat("Smoothing parameter lambda     :", format(signif(x$lambda)), "\n") 
}
