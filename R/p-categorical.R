
#dyn.load("/Users/stasinom/Dropbox/gamlss/R-code/Penalised-Smoothers/GEND/genD.so")
#is.loaded("genD")
# function for merging levels of a factors 
# Based on ideas of Gerhard Tutz (IWSM 2014)
# Paul Eilers, Mikis Stasinopoulos and Marco Enea
#-------------------------------------------------------------------------------
# QUESTIONS
# 1) Is the method of estimating lqmbda correct?
#-------------------------------------------------------------------------------
# TO DO
#  to determine kappa Probably OK
#  Adjust y or adjust kappa? probrly OK
#  fixing df's is not working for Lp not equal to 2 is fixed but slow 
#-------------------------------------------------------------------------------
# this is a new experimental function to fit categorical variables
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#===============================================================================
################################################################################
#===============================================================================
#==============================================================================
################################################################################
#===============================================================================
pcat<- function(fac, df = NULL, 
                lambda = NULL,
                method = c("ML","GAIC"),
                 start = .001,  
                    Lp = 0,
                 kappa = 1e-5, 
                  iter = 100,  # no of iterations 
                c.crit = 1.0e-4,
                     k = 2) 
{
     scall <- deparse(sys.call(), width.cutoff = 500L)
    method <- match.arg(method)
  # check for standarized matric
        X <- model.matrix(~fac - 1)
        p <- ncol(X)
        n <- nrow(X)
# Penalty matrix for all possible differences
#-------------------------------------------------------------------------------
# Marco's new function to generate D
genD <- function(nc){
      nr <- nc*(nc-1)/2
  zero_vect <- rep(0,nc*nr)
       vect <- .C("genD",as.integer(nc),as.integer(zero_vect), PACKAGE="gamlss")[[2]]
          D <- matrix(vect,nr,nc,byrow=TRUE)
#  rownames(D) <- rep("d",nr) #Do you really need it??
  D
}
#-------------------------------------------------------------------------------
        D <- genD(p)
# the ald function
#    for (i in 2:p) 
#       {
#        for (j in 1:(i - 1)) 
#        {
#      #  	cat(i,j, "\n")
#         d <- rep(0, p)
#      d[i] <- 1
#      d[j] <- -1
#         D <- rbind(D, d)
#        }
#       }
      nd <- nrow(D)
# finish penalty
  if(!is.null(lambda))
  {
    if(lambda<0)
    {
      lambda <- 1
      warning(paste("lambda was negative; have used ", lambda))
    }
  } 
  if(!is.null(df))
  {
    if (df > p)  warning("the df are set to ", p, "\n") 
    df <- if (df > p)  p  else  df  
    if (df < 0)  warning("the df are set to 0") 
    df <- if (df < 0)  0  else  df  
    if (df==p) lambda <- 0  
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
  sl <- sample(letters, 4)
  fourLetters <- paste(paste(paste(sl[1], sl[2], sep=""), sl[3], sep=""),sl[4], sep="")
  startLambdaName <- paste("start.Lambda",fourLetters, sep=".")
  ## put the starting values in the gamlss()environment
  #--------
  assign(startLambdaName, start, envir=gamlss.environment)
  #--------
                         x <- rep(0, n)
              attr(x, "X") <- X
           attr(x, "call") <- substitute(gamlss.pcat(data[[scall]], z, w)) 
              attr(x, "D") <- D
         attr(x, "lambda") <- lambda
         attr(x, "df")     <- df
         attr(x, "method") <- method
            attr(x, "Lp")  <- Lp
         attr(x, "kappa")  <- kappa
         attr(x, "iter")   <- iter
         attr(x, "c.crit") <- c.crit
    #     attr(x,"order")  <- order
         attr(x, "start")  <- start 
             attr(x, "k")  <- k 
    attr(x, "gamlss.env")  <- gamlss.environment
  attr(x, "NameForLambda") <- startLambdaName
                  class(x) <- c("smooth", class(x))  
  x
}
#-------------------------------------------------------------------------------
gamlss.pcat <- function(x, y, w, xeval = NULL, ...)
{  
#-------------------------------------------------------------------------------
# local functions
#-------------------------------------------------------------------------------
  regpen <- function(lambda)
  {
    for (it in 1: iter) 
    {  
         RD <- rbind(R,sqrt(lambda)*sqrt(omega.)*D,P0) #
      svdRD <- svd(RD)  
       rank <- sum(svdRD$d>max(svdRD$d)*.Machine$double.eps^.8)
         U1 <- svdRD$u[1:p,1:rank]  
         y1 <- t(U1)%*%Qy
       beta <- svdRD$v[,1:rank] %*%(y1/svdRD$d[1:rank])
         dm <- max(abs(sm - beta))
         sm <- beta
          u <- D %*% beta # should I multiply sqrt(omega.)*
     omega. <- c(1 / (abs(u)^(2-Lp) + kappa ^ 2)) # L_p
     assign("Omega", value=omega., pos = -1)
  #   Omega <- diag(omega.)
     if (dm < c.crit) break # if difference small stop    
    }
        HH <- (svdRD$u)[1:p,1:rank]%*%t(svdRD$u[1:p,1:rank])
       edf <- sum(diag(HH))
        fv <- X%*%beta
       out <-  list(fv=fv, beta=beta, edf=edf, omega=omega.)  
  }
#------------------------------------------------------------------
#------------------------------------------------------------------
# ## function to find lambdas miimizing the local GAIC        
fnGAIC <- function(lambda, k)
{
  fit <- regpen(lambda)
  fv <- fit$fv         
  GAIC <- sum(w*(y-fv)^2)+k*fit$edf 
  # cat("GAIC", GAIC, "\n")
  GAIC   
}
#------------------------------------------------------------------
#-----------------------------------------------------------------
# ## function to find the lambdas which minimise the local GCV 
# fnGCV <- function(lambda, k)
# {
#   I.lambda.D <- (1+lambda*UDU$values)
#   edf <- sum(1/I.lambda.D)
#   y_Hy2 <- y.y-2*sum((yy^2)/I.lambda.D)+sum((yy^2)/((I.lambda.D)^2))
#   GCV <- (n*y_Hy2)/(n-k*edf)^2
#   GCV
# }  
# 

#-----------------------------------------------------------------
#-----------------------------------------------------------------
# main function starts here 
# ----------------------------------------------------------------
     X <-  if (is.null(xeval)) as.matrix(attr(x,"X"))
           else as.matrix(attr(x,"X"))[seq(1,length(y)),]
              D <- as.matrix(attr(x,"D"))
         lambda <- as.vector(attr(x,"lambda"))
             df <- as.vector(attr(x,"df"))  
             Lp <- as.vector(attr(x,"Lp")) 
          kappa <- as.vector(attr(x,"kappa")) 
           iter <- as.vector(attr(x,"iter")) 
              k <- as.vector(attr(x,"k")) 
         c.crit <- as.vector(attr(x,"c.crit"))
         method <- as.character(attr(x,"method")) 
     gamlss.env <- as.environment(attr(x, "gamlss.env"))
startLambdaName <- as.character(attr(x, "NameForLambda")) 
              N <- sum(w!=0) # DS+FDB 3-2-14
              n <- nrow(X)
              p <- ncol(X)
             aN <- nrow(D)
             ap <- ncol(D)  
            qrX <- qr(sqrt(w)*X, tol=.Machine$double.eps^.8)  
              R <- qr.R(qrX)
              Q <- qr.Q(qrX) 
             Qy <- t(Q)%*%(sqrt(w)*y)
  # 
  if(p!=ap) stop("the dimensions of the penalty matrix and of the design matrix are incompatible")
            P0 <- diag(p) * 1e-6
  ## starting values for smoothing function 
            sm <- rep(0, p) 
        omega. <- rep(1, aN)
          tau2 <- sig2 <- NULL
# now the action depends on the values of lambda and df
#--------------------------------------------------------------
   lambdaS <- get(startLambdaName, envir=gamlss.env) ## geting the starting value
if (lambdaS>=1e+07) lambda <- 1e+07 # MS 19-4-12
if (lambdaS<=1e-07) lambda <- 1e-07 # MS 19-4-12
# cat(lambda, "\n")
# case 1: if lambda is known just fit --------------------------------------------
if (is.null(df)&&!is.null(lambda)||!is.null(df)&&!is.null(lambda))
{
  fit <- regpen(lambda)
   fv <- fit$fv        
} # case 2: if lambda is estimated -------------------------------------------- 
  else if (is.null(df)&&is.null(lambda)) 
  {
    lambda <- lambdaS  # MS 19-4-12
    switch(method,
"ML"={
    for (it in 1:20) 
    {
        fit  <- regpen(lambda)
      gamma. <-  (sqrt(fit$omega)*D) %*% as.vector(fit$beta) # MS ?????
    #  plot(gamma.)
          fv <- X %*% fit$beta
        sig2 <- sum(w * (y - fv) ^ 2) / (N - fit$edf)
        tau2 <- sum(gamma. ^ 2) / (fit$edf)# Tuesday, March 17, 2009 at 11:57
      if(tau2<1e-7) tau2 <- 1.0e-7 # MS 19-4-12
  lambda.old <- lambda
      lambda <- sig2 / tau2
      if (abs(lambda-lambda.old) < 0.0001||lambda>100000) break
   #   cat("lambda sig2e sig2b",lambda,sig2,tau2, '\n')
    }
  },
"GAIC"=  #--------------------------------------------------------------- GAIC
{
  lambda <- nlminb(lambda, fnGAIC,  lower = 1.0e-7, upper = 1.0e7, k=k)$par 
  
    fit  <- regpen(lambda)
      fv <- fit$fv     
  assign(startLambdaName, lambda, envir=gamlss.env)
 }#,
# "GCV"={   #-------------------------------------------------------------- GCV
#   # 
#   wy <- sqrt(w)*y
#   y.y <- sum(wy^2)
#   Rinv <- solve(R)
#   S <- t(D)%*%D
#   UDU <- eigen(t(Rinv)%*%S%*%Rinv)
#   yy <- t(UDU$vectors)%*%Qy #t(qr.Q(QR))%*%wy
#   lambda <- nlminb(lambda, fnGCV,  lower = 1.0e-7, upper = 1.0e7, k=k)$par
#   fit <- regpen(lambda)
#   fv <- X %*% fit$beta     
#   assign(startLambdaName, lambda, envir=gamlss.env) 
# }
)
  }
  else # case 3 : if df are required
  { 
    #method 2 from Simon Wood (2006) pages 210-211, and 360 
    ## local function to get df using eigen values
    edf1_df <- function(lambda)
    {
      if (Lp!=2)
      {
        fit  <- regpen(lambda)
         out <- (fit$edf-df)
      } else
      {
        edf <-  sum(1/(1+lambda*UDU$values))
        out <- (edf-df)
      }
    out
    }   
    Rinv <- solve(R)
     S   <- t(D)%*%D
     UDU <- eigen(t(Rinv)%*%S%*%Rinv)           
  lambda <- if (sign(edf1_df(0))==sign(edf1_df(100000))) 100000  # in case they have the some sign
    else  uniroot(edf1_df, c(0,100000))$root
    # if (any(class(lambda)%in%"try-error")) {lambda<-100000}  
   fit  <- regpen(lambda)
     fv <- fit$fv
  }      
  waug <- as.vector(c(w, rep(1,nrow(D))))
  xaug <- as.matrix(rbind(X,sqrt(lambda)*(sqrt(fit$omega)*D))) # MS ?????
   lev <- hat(sqrt(waug)*xaug,intercept=FALSE)[1:n] # get the hat matrix
 # lev <- (lev-.hat.WX(w,rep(1,n)))
  #lev <- (lev-.hat.WX(w,x)) # substract  the linear since is allready fitted 
  var <- lev/w # the variance of the smoother
coefSmo <- list(coef = fit$beta, 
              lambda = lambda, 
                 edf = fit$edf, 
               sigb2 = tau2, 
               sige2 = sig2,
                  fv = fv,  # MS ????
             factor  = factor(zapsmall(ifelse(abs(fv)<0.0000001,0,fv), digits=3)), # ????
                  se = sqrt(var),
                  Lp = Lp)
class(coefSmo) <- "pcat"
  #----------------------------------------------
  if (is.null(xeval))
  {
    list(fitted.values=fv, residuals=y-fv, var=var, nl.df =fit$edf-1,
         lambda=lambda, coefSmo=coefSmo )
  }
  else 
  {
    ll <- dim(as.matrix(attr(x,"X")))[1]
    nx <- as.matrix(attr(x,"X"))[seq(length(y)+1,ll),]
    pred <- drop(nx %*% fit$beta)
    pred
  }    
}
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
plot.pcat <- function(x, ...)
{
  plot(x$coef, type="h", xlab="levels", ylab="coefficients", 
       main=paste("Lp =", paste(x$Lp), sep=" "))
  abline(h=0)
}
#-------------------------------------------------------------------------------
coef.pcat <- function(object, ...)
{
as.vector(object$coef)
}
#-------------------------------------------------------------------------------
print.pcat <- function (x, digits = max(3, getOption("digits") - 3), ...) 
{   
  cat("Randon effects fit using the gamlss function pcat() \n")
  cat("Degrees of Freedom for the fit :", x$edf, "\n")
  cat("Random effect parameter sigma_b:", format(signif(x$sigb)), "\n")  
  cat("Smoothing parameter lambda     :", format(signif(x$lambda)), "\n") 
}
#-------------------------------------------------------------------------------
fitted.pcat<- function(object,...)
{
  as.vector(object$fv)
}
#===============================================================================
################################################################################
#===============================================================================
plotDF <- function(y, 
                        factor = NULL, 
                        formula = NULL, 
                        data, 
                        along=seq(0,nlevels(factor)), 
                        kappa=0.000001, 
                        Lp=0, 
                        ...)
{
  #if (is.null(factor)) stop("the factor argument is needed")
  if (is.null(formula))
  {
    ylab <- deparse(substitute(y))
    xlab <- deparse(substitute(factor))
    y <- if (!is.null(data)) get(deparse(substitute(y)), envir=as.environment(data)) else y
    factor <- if (!is.null(data)) get(deparse(substitute(factor)), envir=as.environment(data)) else factor 
    kl <- length(along)
    matdffit <-  matrix(NA, nrow =  kl, ncol =  nlevels(factor))
    for (i in 1:kl)
    {
      cat("df ", along[i], "\n")
      mm <- gamlss(y~pcat(factor, df=along[i], kappa=kappa, Lp=Lp),  trace=F, gd.tol=Inf, ...)
      matdffit[i,] <- coef(getSmo(mm))
    }
  } else 
  {

    if (is.null(data)) stop("the data argument is needed with formula")
      xlab <- deparse(substitute(factor))
    factor <- get(deparse(substitute(factor)), envir=as.environment(data))  
        kl <- with(data, length(along)) 
    matdffit <-  matrix(NA, nrow =  kl, ncol =  nlevels(factor))
   form <- as.formula( paste(paste(paste(formula[2],formula[1]), formula[3]),  paste(paste("pcat(",xlab, sep=""), ", df=along[i], kappa=kappa, Lp=Lp)"), sep="+"))#, env=as.environment(data)
   term <-paste(paste("pcat(",xlab, sep=""),", df=along[i], kappa=kappa, Lp = Lp)", sep="")
   for (i in 1:kl)
   {
     cat("df ", along[i], "\n")
     mm <-  gamlss(form,  data=data,trace=F, gd.tol=Inf, ...)
     matdffit[i,] <-  coef(getSmo(mm, which=dim(mm$mu.s)[2]))
   }
  
  }  
  matplot( along, matdffit, type="l", xlim=c(0,nlevels(factor)+1), xlab="df")
  text(rep(along[kl]+0.5,  nlevels(factor)+0.1), matdffit[kl, ], levels(factor))
}
#-------------------------------------------------------------------------------
plotLambda  <- function(y, 
                        factor = NULL, 
                        formula = NULL, 
                        data, 
                        along= seq(-2,2, .1), 
                        kappa=0.000001, 
                        Lp=0, 
                        ...)
{
  #if (is.null(factor)) stop("the factor argument is needed")
  if (is.null(formula)) # if not formula 
  {
    ylab <- deparse(substitute(y))
    xlab <- deparse(substitute(factor))
    y <- if (!is.null(data)) get(deparse(substitute(y)), envir=as.environment(data)) else y
    factor <- if (!is.null(data)) get(deparse(substitute(factor)), envir=as.environment(data)) else factor 
    kl <- length(along)
    matdffit <-  matrix(NA, nrow =  kl, ncol =  nlevels(factor))
    for (i in 1:kl)
    {
      cat("lambda ", along[i], "\n")
      mm <- gamlss(y~pcat(factor, lambda=exp(along[i]), kappa=kappa, Lp=Lp),  trace=F, gd.tol=Inf, ...)
      matdffit[i,] <- coef(getSmo(mm))
    }
  } else  # if formula 
  {
    if (is.null(data)) stop("the data argument is needed with formula")
      xlab <- deparse(substitute(factor))
    factor <- get(deparse(substitute(factor)), envir=as.environment(data))  
        kl <- with(data, length(along)) 
    matdffit <-  matrix(NA, nrow =  kl, ncol =  nlevels(factor))
    form <- as.formula( paste(paste(paste(formula[2],formula[1]), formula[3]),  paste(paste("pcat(",xlab, sep=""), ", lambda=exp(along[i]), kappa=kappa, Lp=Lp)"), sep="+"))#, env=as.environment(data)
#    term <-paste(paste("pcat(",xlab, sep=""),", lambda=exp(along[i]), kappa=kappa, Lp = Lp)", sep="")
    for (i in 1:kl)
    {
      cat("lambda ", along[i], "\n")
      mm <-  gamlss(form,  data=data,trace=F, gd.tol=Inf, ...)
      matdffit[i,] <-  coef(getSmo(mm, which=dim(mm$mu.s)[2]))
    }
  }  
  matplot(along, matdffit, type="l", xlab="log-lambda")
  text(rep(along[1]-0.1,  nlevels(factor)+0.1), matdffit[1, ], levels(factor))
}

#-------------------------------------------------------------------------------
# genD <-  function(nl)
# {
#   D <- NULL
#   for (i in 2:nl) 
#   {
#     for (j in 1:(i - 1)) 
#     {
#       cat(i,j, "\n")
#       d <- rep(0, nl)
#       d[i] <- 1
#       d[j] <- -1
#       D <- rbind(D, d)
#     }
#   } 
#   D
# }
