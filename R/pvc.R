## this is part of the implementation of the Penalized B-splines smoothers
## Paul Eilers, Mikis Stasinopoulos and Bob Rigby
## last modified Saturday, May 8, 2010 
## the pvc() function is a varying coefficients function
#----------------------------------------------------------------------------------------
pvc<-function(x, df = NULL, lambda = NULL, by = NULL, control=pvc.control(...), ...) 
{
## this function is based on Paul Eilers' penalised beta regression splines function
## lambda : is the smoothing parameter
## df : are the effective df's
## if both lambda=NULL  and df=NULL then lambda is estimated using the different method
## methods are "ML", "ML-1",  "GAIC" and "GCV"  
## if df is set to number but lambda is NULL then df are used for smoothing
## if lambda is set to a number (whether df=NULL  or not) lambda is used for smoothing
## the knots are defined in equal space unless quantiles=TRUE is set in the control function
# ---------------------------------------------------
## local function
## creates the basis for p-splines
## Paul Eilers' function
#----------------------------------------------------
 bbase <- function(x, xl, xr, ndx, deg, quantiles=FALSE)
  {
 tpower <- function(x, t, p)
# Truncated p-th power function
    (x - t) ^ p * (x > t)
# DS xl= min, xr=max, ndx= number of points within 
# Construct B-spline basis
     dx <- (xr - xl) / ndx # DS increment 
  knots <- if (quantiles==TRUE)
               { c(seq(xl-deg*dx, xl, dx),quantile(x, prob=seq(0, 1, length=ndx)),seq(xr, xr+deg*dx, dx))}
          else {seq(xl - deg * dx, xr + deg * dx, by = dx)} 
      P <- outer(x, knots, tpower, deg)# calculate the power in the knots
      n <- dim(P)[2]
      D <- diff(diag(n), diff = deg + 1) / (gamma(deg + 1) * dx ^ deg) # 
      B <- (-1) ^ (deg + 1) * P %*% t(D) 
      B 
  }
#---------------------------------------------------
# the main function starts here
         scall <- deparse(sys.call())
            lx <- length(x)
 control$inter <- if (lx<100) min(c(10,control$inter)) else control$inter # this is to prevent singularities when length(x) is small
            xl <- min(x)
            xr <- max(x)
          xmax <- xr + 0.01 * (xr - xl)
          xmin <- xl - 0.01 * (xr - xl)   
             X <- bbase(x, xmin, xmax, control$inter, control$degree, control$quantiles) # create the basis
             r <- ncol(X)
             D <- if(control$order==0) diag(r) else diff(diag(r), diff=control$order) # the penalty matrix
## here we get the gamlss environment and a random name to save
## the starting values for lambda within gamlss()
## get gamlss environment
#--------
            rexpr <- regexpr("gamlss",sys.calls())
   for (i in length(rexpr):1)
   { 
    position <- i 
    if (rexpr[i]==1) break
   }
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
if(is.null(by))
   {
     by.var <- NULL
       xvar <- x   
     if(!is.null(df)) # degrees of freedom
             {
             if (df[1]>(dim(X)[2]-2)) 
              {df <- 3;  
              warning("The df's exceed the number of columns of the design matrix", "\n",  "   they are set to 3") }
              df <- if (any(df < 1))  1  else  df+2
              if (any(df < 1) ) warning("the df are set to 1")    
             }
   }
   else
   {
     if (is(by, "factor")) 
      { by.var <- by 
      starting <- rep(control$start, nlevels(by))
      assign(startLambdaName, starting, envir=gamlss.environment)
       if(!is.null(df[1])) # degrees of freedom
             {
             if (length(df)!=nlevels(by)) df <- rep(df[1], nlevels(by))
             if (any(df>(dim(X)[2]-2))) 
              {df <- ifelse((df>(dim(X)[2]-2)),3,df)  
              warning("The df's exceed the number of columns of the design matrix", "\n",  "   they are set to 3") }
              df <- ifelse( df < 1,  1 ,  df+2)
              if (any(df < 1))  warning("the df are set to 1")    
             }
            }
     else
      { by.var <- by-mean(by)
 if(!is.null(df)) # degrees of freedom
             {
             if (df>(dim(X)[2]-2)) 
              {df <- 3;  
              warning("The df's exceed the number of columns of the design matrix", "\n",  "   they are set to 3") }
              df <- if (df[1] < 1)  1  else  df+2
              if (df[1] < 1)  warning("the df are set to 1")    
             }      	
      }
    xvar <- model.matrix(~by.var*x, contrast="")  
   }
      attr(xvar, "control")       <- control
      attr(xvar, "D")             <- D
      attr(xvar, "X")             <- X
      attr(xvar, "df")            <- df 
      attr(xvar, "by")            <- by.var
      attr(xvar, "call")          <- substitute(gamlss.pvc(data[[scall]], z, w)) 
      attr(xvar, "lambda")        <- lambda
      attr(xvar, "gamlss.env")    <- gamlss.environment
      attr(xvar, "NameForLambda") <- startLambdaName
      attr(xvar, "class")         <- "smooth"
      xvar
}
#----------------------------------------------------------------------------------------
# control function for pvc()
##---------------------------------------------------------------------------------------
pvc.control <- function(inter = 20, degree= 3, order = 2, start=10, quantiles=FALSE, 
                       method=c("ML","GAIC", "GCV", "EM", "ML-1"), k=2, ...)
{ 
##  Control function for pb()
##  MS  Tuesday, March 24, 2009
## inter : is the number of equal space intervals in x 
## (unless quantiles = TRUE is used in which case the points will be at the quantiles values of x) 
## degree: is the degree of the polynomial 
## order refers to differences in the penalty for the coeficients 
## order = 0 : white noise random effects
## order = 1 : random walk
## order = 2 : random walk of order 2
## order = 3 : random walk of order 3
# inter = 20, degree= 3, order = 2, start=10, quantiles=FALSE, method="loML"
        if(inter <= 0) {
warning("the value of inter supplied is less than 0, the value of 10 was used instead")
                inter <- 10 }
        if(degree <= 0) {
warning("the value of degree supplied is less than zero or negative the default value of 3 was used instead")
                degree <- 3}                
        if(order < 0) {
warning("the value of order supplied is zero or negative the default value of 2 was used instead")
                order <- 2}
        if(k <= 0) {
warning("the value of GAIC/GCV penalty supplied is less than zero the default value of 2 was used instead")
                k <- 2}   
method <- match.arg(method)                          
        list(inter = inter, degree = degree,  order = order, start=start, 
                   quantiles = as.logical(quantiles)[1], method= method, k=k)
}
#----------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------
gamlss.pvc <- function(x, y, w, xeval = NULL, ...)
{
# -------------------------------------------------- 
## functions within
## local function for simple penalized regression
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
#--------------------------------------------------
## local function for penalised regression covariance model ~x+factor
## using expanding data and lm.fit
regpenByFactor <- function(y, X, w, fac, lambda, D) 
 {
          BX <- FbyX(fac, X)   # get the interaction design matrix 
          BD <- ExpandD(fac, D, lambda) # get the expanded penalty matrix   
        waug <- as.vector(c(w, rep(1,nlev*p2)))
        yaug <- as.vector(c(y, rep(0,nlev*p2)))
        xaug <- as.matrix(rbind(BX,BD))
        fit1 <- lm.wfit(xaug,yaug,w=waug,method="qr") 
         lev <- hat(sqrt(waug)*xaug,intercept=FALSE)[1:n] # get the hat matrix
         fit <- list(fv=fitted(fit1)[1:n], coef=coef(fit1), edf=sum(lev))
        return(fit)  
 }
#--------------------------------------------------
## local function for penalised regression covariance model ~x+factor
## using matrix operators
regpenByFactor1 <- function(y, X, w, fac, lambda, D) # 
 {
          BX <- FbyX(fac, X)   # get the interaction design matrix 
          BD <- ExpandD(fac, D, lambda) # get the expanded penalty matrix   
           G <- t(BD) %*% BD
          XW <- w * BX
         XWX <- t(XW) %*% BX
        beta <- solve(XWX + G, t(XW) %*% y)
          fv <- BX %*%beta
           H <- solve(XWX + G, XWX)
         fit <- list(fv=fv, beta = beta, edf = sum(diag(H)))
        return(fit)  
 }    
#--------------------------------------------------
## similar to regpen  but with saving extra quantities
regpenEM <- function(y, X, w, lambda, order, D)
  {
             G <- lambda * t(D) %*% D
            XW <- w * X
           XWX <- t(XW) %*% X
          beta <- solve(XWX + G, t(XW) %*% y)
            fv <- X %*%beta
             H <- solve(XWX + G, XWX)
             V <- solve(XWX + G)
           fit <- list(beta = beta, edf = sum(diag(H)), V=V)
  return(fit)  
  }
#--------------------------------------------------
## local function to find lambda minimizing the local GAIC
     fnGAIC <- function(lambda, k)
    {
       fit <- regpen(y=y, X=X, w=w, lambda=lambda, D)
        fv <- X %*% fit$beta         
      GAIC <- sum(w*(y-fv)^2)+k*fit$edf 
    # cat("GAIC", GAIC, "\n")
      GAIC   
    }
#
#--------------------------------------------------
## local function to find lambdas miimizing the local GAIC
## here is for model ~X+factor
    fnGAICforFactor <- function(lambda, k)
    {
       fit <- regpenByFactor1(y=yvar, X=X, w=w, fac=by.var, lambda=lambda, D=D) 
        fv <- fit$fv         
      GAIC <- sum(w*(yvar-fv)^2)+k*fit$edf 
    # cat("GAIC", GAIC, "\n")
      GAIC   
    }
#
#--------------------------------------------------
## local function to find lambda minimising the local GCV 
      fnGCV <- function(lambda, k)
           {
    I.lambda.D <- (1+lambda*UDU$values)
           edf <- sum(1/I.lambda.D)
         y_Hy2 <- y.y-2*sum((yy^2)/I.lambda.D)+sum((yy^2)/((I.lambda.D)^2))
           GCV <- (n*y_Hy2)/(n-k*edf)
           GCV
           }  
#--------------------------------------------------
## local function to find lambdas minimising the local GCV
## this is the case where the model is ~x+factor  
      fnGCVforFactor <- function(lambda, k)
           {
           fac <- as.factor(colnames(S))
          lambdas <- model.matrix(~fac-1)%*%lambda  
    I.lambda.D <- (1+lambdas*UDU$values)
           edf <- sum(1/I.lambda.D)
         y_Hy2 <- y.y-2*sum((yy^2)/I.lambda.D)+sum((yy^2)/((I.lambda.D)^2))
           GCV <- (n*y_Hy2)/(n-k*edf)
           GCV
           }  
#--------------------------------------------------
#--------------------------------------------------
## local function to get df using eigen values
    edf1_df <- function(lambda)
           {
           edf <-  sum(1/(1+lambda*UDU$values))
           (edf-df)
           }
#-------------------------------------------------
## local function to  create the interaction desing matrix of the kind |X1 0|
         FbyX<- function(fac, X)  #                                |0 X2|                        
           { # this function is to create a factor * matrix interaction design matrix      
               A  <- model.matrix(~fac-1, contrast="") #get the dummy matrix of the factor
             nlev <- nlevels(fac) 
              for (i in 1:nlev)
                 {
                  NX <- A[,i]* X  # interaction with each column                 
                  XX <- if (i==1) NX else cbind(XX,NX) # put it together
                 }
                XX 
            }
#------------------------------------------------- 
## this function creates a general penalty function
## set lambda=1 if you just want the D's 
        ExpandD <- function(fac,D, lambda)    #  i.e.  | sqrt(lambda1)D      0         |
        { #creating a large penalty matrix             |          0     sqrt(lambda2)D |
          pos1 <- 1
          pos2 <- p2
          pos3 <- 1
          pos4 <- p
            BD <- matrix(0, nrow=nlev*p2, ncol=(nlev*p))
         cname <- rep(0, nlev*p)
           for (j in 1:nlev)
             {
              prow <- (pos1:pos2)
              pcol <- (pos3:pos4)
              pos1 <- pos2+1
              pos2 <- pos1+p2-1
              pos3 <- pos4+1
              pos4 <- pos3+p-1
                DD <- sqrt(lambda[j])*D
    BD[prow, pcol] <- DD
       cname[pcol] <- j 
              }
      colnames(BD) <- cname
         BD
          }
#------------------------------------------------------------------
#------------------------------------------------------------------
#------------------------------------------------------------------
# the main function starts here
# get the attributes
 is.Factor <- FALSE # whether by is factor or variate
         X <-  if (is.null(xeval)) as.matrix(attr(x,"X")) #the trick is for prediction
              else  as.matrix(attr(x,"X"))[seq(1,length(y)),]
         D <- as.matrix(attr(x,"D")) # penalty
    lambda <- as.vector(attr(x,"lambda")) # lambda
        df <- as.vector(attr(x,"df")) # degrees of freedom
    by.var <- if (is.null(xeval)) attr(x,"by")
              else attr(x,"by")[seq(1,length(y))] 
            if (!is.null(by.var)) 
                 {
                 if (is.factor(by.var))
                    {  
                     by.var <- as.factor(by.var)
                  is.Factor <- TRUE
                    }
                 else
                    {  
                     by.var <- as.vector(by.var)
                    }
                 } 
   control <- as.list(attr(x, "control")) 
gamlss.env <- as.environment(attr(x, "gamlss.env"))
startLambdaName <- as.character(attr(x, "NameForLambda")) 
     order <- control$order # the order of the penalty matrix
         n <- nrow(X) # the no of observations
         p <- ncol(D) # the rows of the penalty matrix
        p2 <- nrow(D)
     tau2  <- sig2 <- NULL
        w  <- if (is.null(by.var)) w
              else 
               {
                if (is.factor(by.var)) w
                else w*(by.var^2)
                }
     yvar  <-  if (is.null(by.var)) y
              else 
               {
                if (is.factor(by.var)) y
                else (y/ifelse(by.var==0,0.0001,by.var)) 
                }
## we need to know whether by.var is a factor or not
## if is not a factor fit the models with weights --------------------VARIABLE-VARIABLE    
if (!is.Factor)
 {
## now the action depends on the values of lambda and df
##--------------------------------------------------------------------  
## case 1: if lambda is known just fit
  if (is.null(df)&&!is.null(lambda)||!is.null(df)&&!is.null(lambda))
  {
           fit <- regpen(yvar, X, w, lambda,  D)
            fv <-  X %*% fit$beta          
  } # case 2: if lambda is estimated ----------------------------------
  else if (is.null(df)&&is.null(lambda)) # estimate lambda 
  { #   
   # cat("----------------------------","\n")
         lambda <-  get(startLambdaName, envir=gamlss.env) ## geting the starting value
   # if ML ----------------------
   switch(control$method,    # which method to use for estimating lambda
   "ML"={  
        for (it in 1:50) 
          {
            fit  <- regpen(yvar, X, w, lambda, D) # fit model
          gamma. <- D %*% as.vector(fit$beta)  # get the gamma differences
              fv <- X %*% fit$beta             # fitted values
            sig2 <- sum(w * (yvar - fv) ^ 2) / (n - fit$edf)
            tau2 <- sum(gamma. ^ 2) / (fit$edf-order)# Monday, March 16, 2009 at 20:00 see LNP page 279
      lambda.old <- lambda
          lambda <- sig2 / tau2 # maybe only 1/tau2 will do since it gives exactly the EM results see LM-1
      if (lambda<1.0e-7) lambda<-1.0e-7 # DS Saturday, April 11, 2009 at 14:18
       #    cat("iter tau2 sig2",it,tau2, sig2, '\n')
      if (abs(lambda-lambda.old) < 1.0e-7||lambda>1.0e7) break
       assign(startLambdaName, lambda, envir=gamlss.env)
      #cat("lambda",lambda, '\n')
          }
        },
   "ML-1"={
        for (it in 1:50) 
          {
            fit  <- regpen(yvar, X, w, lambda, D) # fit model
          gamma. <- D %*% as.vector(fit$beta)  # get the gamma differences
              fv <- X %*% fit$beta             # fitted values
            sig2 <- 1 # sum(w * (y - fv) ^ 2) / (n - fit$edf)
            tau2 <- sum(gamma. ^ 2) / (fit$edf-order)# Monday, March 16, 2009 at 20:00 see LNP page 279
      lambda.old <- lambda
          lambda <- sig2 / tau2 # 1/tau2 
      if (lambda<1.0e-7) lambda<-1.0e-7 # DS Saturday, April 11, 2009
      if (abs(lambda-lambda.old) < 1.0e-7||lambda>1.0e7) break
       assign(startLambdaName, lambda, envir=gamlss.env)
          }
        },
   "EM"={
       for (it in 1:500) 
          {
              fit  <- regpenEM(yvar, X, w, lambda, order, D)
            gamma. <- D %*% as.vector(fit$beta)
            vgamma <- sum(diag(D%*%fit$V%*%t(D))) # this is crucial for estimating the variance of gamma Monday, March 23, 2009
                fv <- X %*% fit$beta
              tau2 <- ((sum(gamma.^ 2))+vgamma)/length(gamma.) 
        lambda.old <- lambda
            lambda <- 1 / tau2
          if (lambda<1.0e-7) lambda<-1.0e-7 # DS Saturday, April 11, 2009    
        #    cat("iter sigma_t^2",it, tau2, "lambda",lambda, '\n')
        if (abs(lambda-lambda.old) < 1.0e-7||lambda>1.0e7) break
          }
     #cat("lambda",lambda, '\n')
       assign(startLambdaName, lambda, envir=gamlss.env)
        },
   "GAIC"=
        {
         lambda <- nlminb(lambda, fnGAIC,  lower = 1.0e-7, upper = 1.0e7, k=control$k)$par 
            fit <- regpen(y=yvar, X=X, w=w, lambda=lambda, D)
             fv <- X %*% fit$beta     
         assign(startLambdaName, lambda, envir=gamlss.env)
        },
   "GCV"={
   # 
            QR <-qr(sqrt(w)*X)
            wy <- sqrt(w)*y
           y.y <- sum(wy^2)
          Rinv <- solve(qr.R(QR))
             S <- t(D)%*%D
           UDU <- eigen(t(Rinv)%*%S%*%Rinv)
            yy <- t(UDU$vectors)%*%t(qr.Q(QR))%*%wy
        lambda <- nlminb(lambda, fnGCV,  lower = 1.0e-7, upper = 1.0e7, k=control$k)$par
           fit <- regpen(y=yvar, X=X, w=w, lambda=lambda, D)
            fv <- X %*% fit$beta     
         assign(startLambdaName, lambda, envir=gamlss.env) 
        })
   }
   else # case 3 : if df are required---------------------------------
   { 
       #method 2 from Simon Wood (2006) pages 210-211, and 360 
            QR <- qr(sqrt(w)*X)
          Rinv <- solve(qr.R(QR))
           S   <- t(D)%*%D
           UDU <- eigen(t(Rinv)%*%S%*%Rinv)           
        lambda <- if (sign(edf1_df(0))==sign(edf1_df(100000))) 100000  # in case they have the some sign
                  else  uniroot(edf1_df, c(0,100000))$root
       # if (any(class(lambda)%in%"try-error")) {lambda<-100000}   
            fit <- regpen(yvar, X, w, lambda, D)
             fv <- X %*% fit$beta
   }#--------------------------------------------------------------------------end of case 3
        waug <- as.vector(c(w, rep(1,nrow(D))))
        xaug <- as.matrix(rbind(X,sqrt(lambda)*D))
         lev <- hat(sqrt(waug)*xaug,intercept=FALSE)[1:n] # get the hat matrix
         lev <- (lev-.hat.WX(w,X)) # subtract  the linear since is already fitted 
         var <- lev/w              # the variance of the smootherz 
   if (is.null(xeval)) # if no prediction 
   {
     if(is.null(by.var))
     {
     list(fitted.values=fv, residuals=y-fv, var=var, nl.df = fit$edf-2, 
           lambda=lambda, coefSmo=list(coef=fit$beta, lambda=lambda, edf=fit$edf, tau2=tau2, sig2=sig2, method=control$method) )
     }
     else
     {
      list(fitted.values=fv*by.var, residuals=y-(fv*by.var), var=var, nl.df=fit$edf-2, 
           lambda=lambda, coefSmo=list(coef=fit$beta, lambda=lambda, edf=fit$edf, tau2=tau2, sig2=sig2, method=control$method) )
     }
   }                            
   else # for prediction 
   {
     if(is.null(by.var)) # if by is null do what you do in pb() 
     {
      ll <- dim(as.matrix(attr(x,"X")))[1]
      nx <- as.matrix(attr(x,"X"))[seq(length(y)+1,ll),]
    pred <- drop(nx %*% fit$beta)
     }
     else # if by is not null you 
     { 
      ll <- dim(as.matrix(attr(x,"X")))[1]
      nx <- as.matrix(attr(x,"X"))[seq(length(y)+1,ll),]
    pred <- drop(nx %*% fit$beta)*by.var
      }
     pred
     }  
 }# if by is NOT a FACTOR ends here -----------------------------------
else  # here is if.Factor==TRUE     -----------------------------------FACTOR---FACTOR
 {#    by  is a FACTOR
           nlev <- nlevels(by.var)
            lev <- levels(by.var)
            fit <- list(nlev)
       coefbeta <- list(nlev)
            fv  <- rep(0, length(y))
           edf  <- rep(0, nlev)
         lambda <- if (length(lambda)!=nlev) rep(lambda, nlev)
                   else lambda
          Index <- 1:length(yvar)
     if (is.null(df[1])&&!is.null(lambda[1])||!is.null(df[1])&&!is.null(lambda[1]))
     { # this is for fix lambdas
     # for (i in 1:nlev)
     #  { 
     #            II <- by.var==lev[i]
     #             In <- Index[II]
     #       fit[[i]] <- regpen(yvar[II], X[II,], w[II], lambda[i],  D)
     #         edf[i] <- fit[[i]]$edf
     #  coefbeta[[i]] <- fit[[i]]$beta 
     #         fv[In] <- X[II,]%*%fit[[i]]$beta      
     #  } 
     ## note that looping or just fitting using regpenByFactor1 is equivalent 
            fit <- regpenByFactor1(yvar, X, w, by.var, lambda, D)   
            fv  <- fit$fv
            edf <- fit$edf 
     } # case 2: if lambda is estimated ----------------------------------
    else if (is.null(df[1])&&is.null(lambda[1])) 
     { #   
      #cat("----------------------------","\n")
         lambda <-  get(startLambdaName, envir=gamlss.env) ## getting the starting value
     # if ML ----------------------
     switch(control$method,
     "ML"={## the method used here is to loop and fit the individual models
           for (i in 1:nlev)
           {
                 II <-  by.var==lev[i]
                 In <- Index[II]
                 #cat("parameter",i,"\n")
                  for (it in 1:50) 
                   {
                fit[[i]] <- regpen(yvar[II], X[II,], w[II], lambda[i],  D)
                  edf[i] <- fit[[i]]$edf
           coefbeta[[i]] <- fit[[i]]$beta 
                     lfv <- X[II,]%*%fit[[i]]$beta              # fitted values
                  fv[In] <- lfv 
                  gamma. <- D %*% as.vector(fit[[i]]$beta)  # get the gamma differences
                    sig2 <- sum(w[II] * (y[II] - lfv) ^ 2) / (sum(II) - fit[[i]]$edf)
                    tau2 <- sum(gamma. ^ 2) / (fit[[i]]$edf-order)# Monday, March 16, 2009 at 20:00 see LNP page 279
              lambda.old <- lambda[i]
               lambda[i] <- sig2 / tau2 # 
               # cat("it lambda", it, lambda, "\n")
              if (lambda[i]<1.0e-7) lambda[i]<-1.0e-7 # DS Saturday, April 11, 2009 at 14:18
              if (lambda[i]>1.0e+7) lambda[i]<-1.0e+7 # DS Wednesday, July 29, 2009 at 19:00
              if (abs(lambda[i]-lambda.old) < 1.0e-7||lambda[i]>1.0e7) break
                    }                  
            }
            assign(startLambdaName, lambda, envir=gamlss.env)       
          },
   "ML-1"={## same as "ML" but with sig2=1
           for (i in 1:nlev)
           {
                 II <-  by.var==lev[i]
                 In <- Index[II]
                  for (it in 1:50) 
                   {
                   # fit  <- regpen(yvar, X, w, lambda, D) # fit model
                fit[[i]] <- regpen(yvar[II], X[II,], w[II], lambda[i],  D)
                  edf[i] <- fit[[i]]$edf
           coefbeta[[i]] <- fit[[i]]$beta 
                  fv[In] <- X[II,]%*%fit[[i]]$beta 
                  gamma. <- D %*% as.vector(fit[[i]]$beta)  # get the gamma differences
                     lfv <- X[II,] %*% fit[[i]]$beta              # fitted values
                    sig2 <- 1 # sum(w[II] * (y[II] - lfv) ^ 2) / (sum(II) - fit[[i]]$edf)
                    tau2 <- sum(gamma. ^ 2) / (fit[[i]]$edf-order)# Monday, March 16, 2009 at 20:00 see LNP page 279
              lambda.old <- lambda[i]
               lambda[i] <- sig2 / tau2 # maybe only 1/tau2 will do since it gives exactly the EM results see LM-1
              if (lambda[i]<1.0e-7) lambda[i]<-1.0e-7 # DS Saturday, April 11, 2009 at 14:18
              if (abs(lambda[i]-lambda.old) < 1.0e-7||lambda[i]>1.0e7) break
                    }
           }
           #     cat("iter tau2 sig2 lambda",it,tau2, sig2, lambda, '\n')
            assign(startLambdaName, lambda, envir=gamlss.env)      
         },
    "EM"={## is not woth persuing this 
          stop("EM method is not implemented try ML-1")
         },
   "GAIC"=
         {  # note that here the expanded model is used so the search for lambda is in mupliple dimension
         lambda <- nlminb(lambda,   fnGAICforFactor,  lower = rep(1.0e-7, nlev), upper = rep(1.0e7, nlev), k=control$k)$par # k=control$k)$par 
            fit <- regpenByFactor1(yvar, X, w, by.var, lambda, D)
             fv <- fit$fv 
            edf <- fit$edf 
       coefbeta <- fit$beta
         assign(startLambdaName, lambda, envir=gamlss.env)
         },
   "GCV"={
   ## note that here the expanded model is used so the search for lambda is in mupliple dimension   
            BX <- FbyX(by.var, X)   # get the interaction design matrix 
            BD <- ExpandD(by.var, D, lambda=rep(1,nlev)) # get the expanded penalty matrix   
            QR <-qr(sqrt(w)*BX)
            wy <- sqrt(w)*y
           y.y <- sum(wy^2)
          Rinv <- solve(qr.R(QR))
             S <- t(BD)%*%BD
           UDU <- eigen(t(Rinv)%*%S%*%Rinv)
            yy <- t(UDU$vectors)%*%t(qr.Q(QR))%*%wy
        lambda <- nlminb(lambda, fnGCVforFactor, lower=rep(1.0e-7,nlev), upper=rep(1.0e7,nlev), k=control$k)$par
           fit <- regpenByFactor1(yvar, X, w, by.var, lambda, D)
            fv <- fit$fv 
           edf <- fit$edf 
      coefbeta <- fit$beta
         assign(startLambdaName, lambda, envir=gamlss.env) 
        })
    }
   else # case 3 : if df are required---------------------------------
   { 
             #--------------------------------------------------
     ## local function to get df using eigen values
     lambda <- if (is.null(lambda))  rep(0, nlev)     
     for (i in 1:nlev)
       {      	
          edf1_df1 <- function(lambda)
           {
           edf <-  sum(1/(1+lambda*UDU$values))
           (edf-df[i])
           }
         #---------------- 
                  II <- by.var==lev[i]
                  In <- Index[II]
                  QR <- qr(sqrt(w[II])*X[II,])
                Rinv <- solve(qr.R(QR))
                 S   <- t(D)%*%D
                 UDU <- eigen(t(Rinv)%*%S%*%Rinv) 
           lambda[i] <- if (sign(edf1_df1(0))==sign(edf1_df1(100000))) 100000  # in case they have the same sign
                       else  uniroot(edf1_df1, c(0,100000))$root
              if (any(class(lambda[i])%in%"try-error")) {lambda[i]<-100000}             
       } 
                fit <- regpenByFactor1(yvar, X, w, by.var, lambda, D)
                 fv <- fit$fv 
                edf <- fit$edf 
           coefbeta <- fit$beta
         assign(startLambdaName, lambda, envir=gamlss.env) 
    }#--------------------------------------------------------------------------end of case 3
     # I need to calculate the hat matrix here for the variance of the smoother
           BX <- FbyX(by.var, X)   # get the interaction design matrix 
          BD <- ExpandD(by.var, D, lambda) # get the expanded oenalty matrix   
        waug <- as.vector(c(w, rep(1,nlev*p2)))
        yaug <- as.vector(c(y, rep(0,nlev*p2)))
        xaug <- as.matrix(rbind(BX,BD))
        #fit1 <- lm.wfit(xaug,yaug,w=waug,method="qr") 
        #fit2 <- regpenByFactor(yvar, X, w, by.var, lambda, D)
        #fit3 <- regpenByFactor1(yvar, X, w, by.var, lambda, D)   
        # the error is here
         lev <- hat(sqrt(waug)*xaug,intercept=FALSE)[1:n] # get the hat matrix
           lev <- (lev-.hat.WX(w,X)) # This has to checked???? 
         var <- lev/w   
     #if (any(is.na(fv)))
     if (is.null(xeval)) # if no prediction 
     {
      list(fitted.values=fv, residuals=y-fv, var=var, nl.df =sum(edf)-(2*nlev),
           lambda=lambda, coefSmo=list(coef=coefbeta, lambda=lambda, edf=sum(edf), EDF=edf, tau2=NULL, sig2=NULL, method=control$method) )
     }
     else
     {# Mikis 27-6-11
        ll <- dim(as.matrix(attr(x,"X")))[1]                # length of original X
        nx <- as.matrix(attr(x,"X"))[seq(length(y)+1,ll),]  # the x-values matrix
       fac <- attr(x,"by")[seq(length(y)+1,ll)]             # get the new values of the factor 
        BX <- FbyX(fac, nx)                                 # get the interaction design matrix 
  longbeta <- as.vector(sapply(fit, function(l) l$beta))    # get the beta as a long vector 
      pred <- drop(BX %*% longbeta)                         # get prediction
      pred
     }       
   }    
}
