#----------------------------------------------------------------------------------------
# this has a new experimental function to fit ridge
# It is not clear how DF's are defined and whether lambda should be calculated differently 
#----------------------------------------------------------------------------------------
# the standarized function
# I used V and R standazation here
#----------------------------------
#standarized <- function(formula, data)
#   {
#        m <- match.call(expand.dots = FALSE)
#   m[[1]] <- as.name("model.frame")
#        m <- eval.parent(m)
#    Terms <- attr(m, "terms")
#        X <- model.matrix(Terms, m)
#        n <- nrow(X)
#        p <- ncol(X)
#    if (Inter <- attr(Terms, "intercept")) 
#    {
#        Xm <- colMeans(X[, -Inter])
#        p <- p - 1
#        X <- X[, -Inter] - rep(Xm, rep(n, p))
#    }
#    else Ym <- Xm <- NA
#    Xscale <- drop(rep(1/n, n) %*% X^2)^0.5
#         X <- X/rep(Xscale, rep(n, p))
#    attr(X, "class") <- c("standarized", "matrix")
#    X
#    }
#----------------------------------------------------------------------------------------
#======================================================================================== 
ri<- function(X, df = NULL, lambda = NULL, order=0, start=10) 
    {
    scall <- deparse(sys.call())
    # check for standarized matric
    if (any(abs(apply(X,2, "mean")>.5))) warning("The design matrix X should be standarized")
   # if (any(abs(apply(X,2, "sd")>.5))) warning("The design matrix X should be standarized")
    p <- ncol(X)
    n <- nrow(X)
    if(!is.null(lambda))
      { 
       if(lambda<0)
          {
           lambda <- 0
           warning(paste("lambda was negative; have used ", lambda))
          }
      } 
       if(!is.null(df))
             {
            df <- if (df < 1)  1  else  df
            if (df < 1)  warning("the df are set to 1")    
             }     
      # this is included here for generality  
       D <- if(order==0) diag(p) else diff(diag(p), diff=order)
       x <- rep(0, n)
    attr(x, "X")      <- X
    attr(x, "call")   <-  substitute(gamlss.ri(data[[scall]], z, w)) 
    attr(x, "D")      <- D
    attr(x, "lambda") <- lambda
    attr(x, "df")     <- df
    attr(x,"order")   <- order
    attr(x, "start")  <- start 
    class(x)          <- c("smooth", class(x))  
    x
   }
#-------------------------------------------------------------------------------------------
gamlss.ri <- function(x, y, w, xeval = NULL, ...)
{  
# a siple penalized regression
regpen <- function(y, X, w, lambda, order, D)
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
#----------------------  
# ---------------------         
        X <-  if (is.null(xeval)) as.matrix(attr(x,"X"))
             else as.matrix(attr(x,"X"))[seq(1,length(y)),]
        D <- as.matrix(attr(x,"D"))
    order <- as.vector(attr(x,"order"))
   lambda <- as.vector(attr(x,"lambda"))
       df <- as.vector(attr(x,"df"))   
        n <- nrow(X)
        p <- ncol(X)
       aN <- nrow(D)
       ap <- ncol(D)   
    if(p!=ap) stop("the dimensions of the augmented matrix and of the design matrix are incompatible")
 if (is.null(df)&&!is.null(lambda)||!is.null(df)&&!is.null(lambda))
 {
          fit <- regpen(y, X, w, lambda, order, D)
           fv <- X %*% fit$beta
           
 } # case 2: if lambda is estimated 
 else if (is.null(df)&&is.null(lambda)) 
 {
        lambda <- as.numeric(attr(x,"start")) 
  for (it in 1:20) 
    {
          fit  <- regpen(y, X, w, lambda, order, D)
   #   nb = ncol(X)
   #   D = diff(diag(nb), diff =   1 )
        gamma. <- D %*% as.vector(fit$beta)
            fv <- X %*% fit$beta
          sig2 <- sum(w * (y - fv) ^ 2) / (n - fit$edf)
          tau2 <- sum(gamma. ^ 2) / (fit$edf-order)# Tuesday, March 17, 2009 at 11:57
    lambda.old <- lambda
        lambda <- sig2 / tau2
    if (abs(lambda-lambda.old) < 0.0001||lambda>100000) break
    #cat("lambda",lambda, '\n')
    }
  }
  else # case 3 : if df are required
  { 
            XW <- w * X
           XWX <- t(XW) %*% X
        edf_df <- function(lambda)
           {
             G <- lambda * t(D) %*% D
             H <- solve(XWX + G, XWX)
           edf <- sum(diag(H))
          # cat("edf", edf, "\n")
           (edf-df)
           }
       lambda <- if (sign(edf_df(0))==sign(edf_df(100000))) 100000
                 else  uniroot(edf_df, c(0,100000))$root
      # if (any(class(lambda)%in%"try-error")) {lambda<-100000}   
           fit <- regpen(y, X, w, lambda, order, D)
            fv <- X %*% fit$beta
  }        
     waug <- as.vector(c(w, rep(1,nrow(D))))
     xaug <- as.matrix(rbind(X,sqrt(lambda)*D))
      lev <- hat(sqrt(waug)*xaug,intercept=FALSE)[1:n] # get the hat matrix
      lev <- (lev-.hat.WX(w,x)) # substract  the linear since is allready fitted 
      var <- lev/w # the variance of the smootherz
     #----------------------------------------------
if (is.null(xeval))
    {
      list(fitted.values=fv, residuals=y-fv, var=var, nl.df =fit$edf-1,
          lambda=lambda, coefSmo=list(coef=fit$beta, lambda=lambda, edf=fit$edf) )
    }
else 
    {
    ll <- dim(as.matrix(attr(x,"X")))[1]
     nx <- as.matrix(attr(x,"X"))[seq(length(y)+1,ll),]
   pred <- drop(nx %*% fit$beta)
   pred
    }    
}
