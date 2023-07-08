## this is the new implementation of the Penalized B-splines smoother
## which allows the smoother to go to the constant rather than linear
## it can be used for model selection
## Mikis Stasinopoulos, Bob Rigby, Vlasisos Voudoris based  Marias Durban idea of  
## using double penalies, one of order 2 and one of order 1.
## The second penalty only applies if the fit has edf close to 1
## it tries to imitate  Simon Woods's model selection in gam() 
## created  Aug-2015 
#-------------------------------------------------------------------------------
pbz <- function(x, df = NULL, lambda = NULL, control=pbz.control(...), ...) 
{
# ------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
## local function
## creates the basis for p-splines
## Paul Eilers' function
#-------------------------------------------------------------------------------
 bbase <- function(x, xl, xr, ndx, deg, quantiles=FALSE)
  {
 tpower <- function(x, t, p)
# Truncated p-th power function
    (x - t) ^ p * (x > t)
# DS xl= min, xr=max, ndx= number of points within 
# Construct B-spline basis
# if quantiles=TRUE use different bases
        dx <- (xr - xl) / ndx # DS increment 
 if (quantiles) # if true use splineDesign
      { 
      knots <-  sort(c(seq(xl-deg*dx, xl, dx),quantile(x, prob=seq(0, 1, length=ndx)), seq(xr, xr+deg*dx, dx))) 
          B <- splineDesign(knots, x = x, outer.ok = TRUE, ord=deg+1)
          return(B)    
      }
     else # if false use Paul's
     { 
      knots <- seq(xl - deg * dx, xr + deg * dx, by = dx)
          P <- outer(x, knots, tpower, deg)# calculate the power in the knots
          n <- dim(P)[2]
          D <- diff(diag(n), diff = deg + 1) / (gamma(deg + 1) * dx ^ deg) # 
          B <- (-1) ^ (deg + 1) * P %*% t(D) 
          B 
     }
  }
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# the main function starts here
            scall <- deparse(sys.call())
      no.dist.val <-  length(table(x))
      if (is.matrix(x)) stop("x is a matric declare it as a vector")
               lx <- length(x)
    control.inter <- if (lx<99) 10 else control$inter # this is to prevent singularities when length(x) is small:change to 99 30-11-11 MS
    control$inter <- if (no.dist.val<=control$inter)  no.dist.val else control.inter 
               xl <- min(x)
               xr <- max(x)
             xmax <- xr + 0.01 * (xr - xl)
             xmin <- xl - 0.01 * (xr - xl)  
##                create the basis
                X <- bbase(x, xmin, xmax, control$inter, control$degree, control$quantiles) # 
                r <- ncol(X)
                D <- diff(diag(r), diff=control$order)
               D1 <- diff(diag(r))              
             if(!is.null(df)) # degrees of freedom
             {
             if (df>(dim(X)[2]-2)) 
              {df <- 3;  
              warning("The df's exceed the number of columns of the design matrix", "\n",  "   they are set to 3") }
              if (df < 0)  warning("the extra df's are set to 0")   
              df <- if (df < 0)  2  else  df+2
             }
##    
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
               sl <- sample(letters, 4)
      fourLetters <- paste(paste(paste(sl[1], sl[2], sep=""), sl[3], sep=""),sl[4], sep="")
  startLambdaName <- paste("start.Lambda",fourLetters, sep=".")
## put the starting values in the gamlss()environment
#--------
   assign(startLambdaName, control$start, envir=gamlss.environment)
#--------
          xvar <- rep(0,length(x)) # zero in the design matrix the rest pass as artributes
      attr(xvar, "control")       <- control
      attr(xvar, "D")             <- D
      attr(xvar, "D1")            <- D1
      attr(xvar, "X")             <- X
      attr(xvar, "x")             <- x
      attr(xvar, "df")            <- df 
      attr(xvar, "call")          <- substitute(gamlss.pbz(data[[scall]], z, w)) 
      attr(xvar, "lambda")        <- lambda
      attr(xvar, "gamlss.env")    <- gamlss.environment
      attr(xvar, "NameForLambda") <- startLambdaName
      attr(xvar, "class")         <- "smooth"
      xvar
}
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# control function for pb()
##------------------------------------------------------------------------------
pbz.control <- function(inter = 20, degree= 3, order = 2, start=c(0.0001,0.0001), 
                      quantiles=FALSE, method=c("ML","GAIC", "GCV"), k=2, 
                      lim=3,  ...)
{ 
        if(inter <= 0) {
warning("the value of inter supplied is less than 0, the value of 10 was used instead")
                inter <- 10 }
        if(degree <= 0) {
warning("the value of degree supplied is less than zero or negative the default value of 3 was used instead")
                degree <- 3}                
        if(order < 2) {
warning("the value of order supplied is less than 2 the default value of 2 was used instead")
                order <- 2}
        if(k <= 0) {
warning("the value of GAIC/GCV penalty supplied is less than zero the default value of 2 was used instead")
                k <- 2}   
method <- match.arg(method)    
        list(inter = inter, degree = degree,  order = order, start=start, 
                   quantiles = as.logical(quantiles)[1], method= method, 
                  k=k, lim=lim)
}
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
gamlss.pbz <- function(x, y, w, xeval = NULL, ...)
{
# ------------------------------------------------------------------------------ 
# functions within
# a simple penalised regression
# this is the original matrix manipulation version but it swiches to QR if it fails
#------------------------------------------------------------------------------
regpen <- function(y, X, w)# original
{
       RD <- rbind(R,sqrt(lambda)*D) # matrix 
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
  #cat("df", df, "\n")
      df1 <- df2 <- 0
 do1order <- FALSE    
  if (df<=lim) # this is crucial for the cut of point
  {
        RD <- rbind(R,sqrt(lambda)*D,sqrt(lambda2)*D1) # matrix  
       RD1 <- rbind(R,sqrt(lambda)*D)
       RD2 <- rbind(R,sqrt(lambda2)*D1)
     svdRD <- svd(RD)                 # U 2pxp D pxp V pxp
    svdRD1 <- svd(RD1)
    svdRD2 <- svd(RD2)
    ##             take only the important values    
      rank <- sum(svdRD$d>max(svdRD$d)*.Machine$double.eps^.8)
        U1 <- svdRD$u[1:p,1:rank]     # U1 p x rank 
    # I am not sure what are consequances in introducing this ???
        y1 <- t(U1)%*%Qy #  t(Q)%*%(sqrt(w)*y)        # rankxp pxn nx1 => rank x 1 vector 
      beta <- svdRD$v[,1:rank] %*%(y1/svdRD$d[1:rank])
        HH <- (svdRD$u)[1:p,1:rank]%*%t(svdRD$u[1:p,1:rank])
       HH1 <- (svdRD1$u)[1:p,1:rank]%*%t(svdRD1$u[1:p,1:rank]) 
       HH2 <- (svdRD2$u)[1:p,1:rank]%*%t(svdRD2$u[1:p,1:rank]) 
        df <- sum(diag(HH))
       df1 <- sum(diag(HH1))
       df2 <- sum(diag(HH2))
  do1order <- TRUE  
  }
  fit <- list(beta = beta, edf = df,  df1=df1, df2=df2, do1order=do1order)
  return(fit)  
}
# #-------------------------------------------------------------------------------
# ## function to find lambdas miimizing the local GAIC        
     fnGAIC <- function(lambda, k)
    {
       fit <- regpen(y=y, X=X, w=w, lambda=lambda)
        fv <- X %*% fit$beta         
      GAIC <- sum(w*(y-fv)^2)+k*fit$edf 
    # cat("GAIC", GAIC, "\n")
      GAIC   
    }
# #-------------------------------------------------------------------------------
# ## function to find the lambdas which minimise the local GCV 
      fnGCV <- function(lambda, k)
           {
    I.lambda.D <- (1+lambda*UDU$values)
           edf <- sum(1/I.lambda.D)
         y_Hy2 <- y.y-2*sum((yy^2)/I.lambda.D)+sum((yy^2)/((I.lambda.D)^2))
           GCV <- (n*y_Hy2)/(n-k*edf)^2
           GCV
           }  
# #-------------------------------------------------------------------------------
    edf1_df <- function(lambda)
           {
           edf <-  sum(1/(1+lambda*UDU$values))
           (edf-df)
           }  
#-------------------------------------------------------------------------------
# the main function starts here
# get the attributes
#w <- ifelse(w>.Machine$double.xmax^.5,.Machine$double.xmax^.5,w )
if (is.null(xeval)) # if no prediction 
    {    
              X <-  if (is.null(xeval)) as.matrix(attr(x,"X")) #the trick is for prediction
                    else  as.matrix(attr(x,"X"))[seq(1,length(y)),]
           xvar <- as.matrix(attr(x,"x")) # x variable
              D <- as.matrix(attr(x,"D")) # main penalty
             D1 <- as.matrix(attr(x,"D1")) # order 1 penalty   
         lambda <- as.vector(attr(x,"lambda")) # lambda 
             df <- as.vector(attr(x,"df")) # degrees of freedom
        control <- as.list(attr(x, "control")) 
     gamlss.env <- as.environment(attr(x, "gamlss.env"))
startLambdaName <- as.character(attr(x, "NameForLambda")) 
          order <- control$order # the order of the penalty matrix
        lambda2 <- control$start[2]
            lim <- control$lim 
              N <- sum(w!=0) # DS+FDB 3-2-14
              n <- nrow(X) # the no of observations
              p <- ncol(D) # the rows of the penalty matrix
            qrX <- qr(sqrt(w)*X, tol=.Machine$double.eps^.8)  
              R <- qr.R(qrX)
              Q <- qr.Q(qrX) 
             Qy <- t(Q)%*%(sqrt(w)*y)
           tau2 <- sig2 <- tau2_2 <- NULL
            df1 <- df2 <- 0
# now the action depends on the values of lambda and df
#------------------------------------------------------------------------------- 
        lambdaS <- get(startLambdaName, envir=gamlss.env) ## geting the starting value
 if (lambdaS[1]>=1e+07) lambda  <- 1e+07 # MS 19-4-12
 if (lambdaS[1]<=1e-07) lambda  <- 1e-07 # MS 19-4-12
 if (lambdaS[2]>=1e+07) lambda2 <- 1e+07 # MS 19-4-12
 if (lambdaS[2]<=1e-07) lambda2  <- 1e-07 # MS 19-4-12
# case 1: if lambda is known just fit -----------------------------------------
 if (is.null(df)&&!is.null(lambda)||!is.null(df)&&!is.null(lambda))
 {
          fit <- regpen(y, X, w)
           fv <- X %*% fit$beta        
 } # case 2: if lambda is estimated -------------------------------------------- 
 else if (is.null(df)&&is.null(lambda)) 
 { #   
  # cat("----------------------------","\n")
     lambda <- lambdaS[1]  
    lambda2 <- lambdaS[2] 
  # if ML --------------------------------------------------------------------ML     
  switch(control$method,
  "ML"={
#    cat("----", "\n")
       for (it in 1:50) 
         {
           fit  <- regpen(y, X, w) # fit model
             fv <- X %*% fit$beta             # fitted values
#cat(fit$do1order, "\n")
          if (fit$do1order)
          {
      gamma. <- D %*% as.vector(fit$beta)  # get the gamma differences
     gamma2. <- D1 %*% as.vector(fit$beta)   
        sig2 <- sum(w * (y - fv) ^ 2) / (N - fit$edf) # DS+FDB 3-2-14
        tau2 <- sum(gamma. ^ 2) / (fit$df1-1)
      tau2_2 <- sum(gamma2. ^ 2) / (fit$df2-1)
                if(tau2<1e-7)     tau2 <- 1.0e-7 # lower limit
                if(tau2_2<1e-7) tau2_2 <- 1.0e-7 
  lambda.old <- lambda
      lambda <- sig2 / tau2 
                if (lambda<1.0e-7) lambda<-1.0e-7 #
                if (lambda>1.0e+7) lambda<-1.0e+7 # DS 29 3 2012
  lambda2.old <- lambda2
     lambda2 <- sig2 / tau2_2     
                if (lambda2<1.0e-7) lambda2<-1.0e-7 #
                if (lambda2>1.0e+7) lambda2<-1.0e+7 # 
     if (abs(lambda-lambda.old) < 1.0e-7||lambda>1.0e10) break     
#     cat("lambda",lambda,lambda2, fit$edf, '\n')
          } else
          { # the standard pb()
      gamma. <- D %*% as.vector(fit$beta)  # get the gamma differences
          fv <- X %*% fit$beta             # fitted values
        sig2 <- sum(w * (y - fv) ^ 2) / (N - fit$edf) # DS+FDB 3-2-14
        tau2 <- sum(gamma. ^ 2) / (fit$edf-order)# see LNP page 279
        if(tau2<1e-7) tau2 <- 1.0e-7 # MS 19-4-12
  lambda.old <- lambda
      lambda <- sig2 / tau2 
                if (lambda<1.0e-7) lambda<-1.0e-7 # 
                if (lambda>1.0e+7) lambda<-1.0e+7 # DS 29 3 2012
  if (abs(lambda-lambda.old) < 1.0e-7||lambda>1.0e10) break         }  
       }
    assign(startLambdaName, c(lambda, lambda2), envir=gamlss.env) 
       },
  "GAIC"=  #--------------------------------------------------------------- GAIC
       {
        lambda <- nlminb(lambda, fnGAIC,  lower = 1.0e-7, upper = 1.0e7, k=control$k)$par 
           fit <- regpen(y=y, X=X, w=w, lambda=lambda)
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
          fit <- regpen(y=y, X=X, w=w, lambda=lambda)
           fv <- X %*% fit$beta     
        assign(startLambdaName, lambda, envir=gamlss.env) 
       })
  }
  else # case 3 : if df are required--------------------------------------------
  { 
         Rinv <- solve(R)
          S   <- t(D)%*%D
          UDU <- eigen(t(Rinv)%*%S%*%Rinv)           
       lambda <- if (sign(edf1_df(0))==sign(edf1_df(100000))) 100000  # in case they have the some sign
                 else  uniroot(edf1_df, c(0,100000))$root
      # if (any(class(lambda)%in%"try-error")) {lambda<-100000}   
           fit <- regpen(y, X, w, lambda)
            fv <- X %*% fit$beta
  }#end of case 3 --------------------------------------------------------------   
#Version 4 -------------------------------------------------- 
        waug <- as.vector(c(w, rep(1,nrow(D))))
        xaug <- as.matrix(rbind(X,sqrt(lambda)*D))
         lev <- hat(sqrt(waug)*xaug,intercept=FALSE)[1:n] # get the hat matrix
#  MIKIS: conclusion is that version 4 the R hat is the faster
#-end -----------------------------------------------------------    
         lev <- (lev-.hat.WX(w,x)) # subtract  the linear since is already fitted 
         var <- lev/w              # the variance of the smoother
         suppressWarnings(Fun <- splinefun(xvar, fv, method="natural"))
coefSmo <- list(   coef = fit$beta,
                     fv = fv, 
                 lambda = lambda, 
                    edf = fit$edf, 
                  sigb2 = tau2, 
                  sige2 = sig2,
                   sigb = if (is.null(tau2)) NA else sqrt(tau2),
                   sige = if (is.null(sig2)) NA else sqrt(sig2),
                 method = control$method,
                    fun = Fun)
class(coefSmo) <- c("pbz", "pb") 
# if (is.null(xeval)) # if no prediction 
#     {
     list(fitted.values=fv, residuals=y-fv, var=var, nl.df =fit$edf-1,
          lambda=lambda, coefSmo=coefSmo )
    }                            
else # for prediction 
    { 
   #   ll <- dim(as.matrix(attr(x,"X")))[1]
   #   nx <- as.matrix(attr(x,"X"))[seq(length(y)+1,ll),]
   # pred <- drop(nx %*% fit$beta) 
      position=0
      rexpr<-regexpr("predict.gamlss",sys.calls())
      for (i in 1:length(rexpr))
        { 
        position <- i 
        if (rexpr[i]==1) break
        }
#cat("New way of prediction in pbz()  (starting from GAMLSS version 5.0-3)", "\n") 
gamlss.environment <- sys.frame(position)
             param <- get("what", envir=gamlss.environment)
            object <- get("object", envir=gamlss.environment)
                TT <- get("TT", envir=gamlss.environment)
         intercept <- get("coef", envir=gamlss.environment)[1]
     smooth.labels <- get("smooth.labels", envir=gamlss.environment)
                ll <- dim(as.matrix(attr(x,"X")))[1]
           newxval <- as.vector(attr(x,"x"))[seq(length(y)+1,ll)]
              pred <- getSmo(object, parameter= param, which=which(smooth.labels==TT))$fun(newxval)
            # pred <- getSmo(object, parameter= param, which=which(TT%in%smooth.labels))$fun(newxval)
            # pred <- getSmo(object, parameter= param, which=which(TT%in%smooth.labels))$fun(xeval)
   pred
    }    
}
#-------------------------------------------------------------------------------
print.pbz  <- function (x, digits = max(3, getOption("digits") - 3), ...) 
{   
  cat("P-spline fit using the gamlss function pbz() \n")
  cat("Degrees of Freedom for the fit :", x$edf, "\n")
  cat("Random effect parameter sigma_b:", format(signif(x$sigb)), "\n")  
  cat("Smoothing parameter lambda     :", format(signif(x$lambda)), "\n") 
}
