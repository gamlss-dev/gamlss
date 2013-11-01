#-------------------------------------------------------------------------------------------
# this has a new experimental function to fit random effects created Thursday, April 10, 2003 at 13:03
# It can do more that the function random() but at a cost in speed
# in this version  expanatory varables can be fitted to the random coefficients
# the arguments 
#    expl= expanatory variables at the first level
#    data1= the data set at the fisrt level
# can be used for that   
# It is not clear how DF's are defined and whether lambda should be calculated differently 
ra <- function(xfactor, xvector = NULL, df = NULL, lambda = NULL, 
               order = 0, estimate = FALSE, expl = NULL, data1 = NULL) 
    {
    scall <- deparse(sys.call())
    # ckecking if the first argument is a factor 
    if(!inherits(xfactor, "factor")) stop("ra() expects a factor as its first argument")
    xfactor <- C(xfactor, rep(0, length(levels(xfactor))), 1) # set one contrasts to zero so it puts zero in the X matrix 
    n.lev <- nlevels(xfactor) 
    if(is.null(df)&is.null(lambda)) stop("Specify at least one: df or lambda")   
    if(!is.null(lambda))
      { 
       if(lambda<0)
          {
           lambda <- 0
           warning(paste("lambda was negative; have used ", lambda))
          }
      }    
     aug <- diag(n.lev)
     if(order != 0) # needs a check whether is ordered factor or not, ra's  make sense only  for order factors
         {
        for(tt in 1:order) 
          {
            aug <- diff(aug)
          }
         } 
    
     if(!is.null(xvector)) 
         {  if(!is.vector(xvector)) stop("ra() expects a vector as its second argument")
                              dm <-model.matrix(~xfactor:xvector-1,contrast="")
         } 
     else   dm<-model.matrix(~xfactor-1,contrast="")
     if(is.null(expl))
         { exp<-NULL  }
     else    
         {
         if(!is(expl,"formula")) stop("the argument expl in ra() expects a formula")
         else exp<-model.matrix(expl,data=data1)
         }
    attr(xfactor, "design.matrix") <- dm
    attr(xfactor, "pen.augment") <- aug
    attr(xfactor, "explanatory") <- exp 
    attr(xfactor, "call") <- substitute(gamlss.ra(data[[scall]], z, w, df = df )) 
    attr(xfactor, "lambda") <- lambda
    attr(xfactor, "df") <- df
    attr(xfactor, "estimate") <- estimate
    attr(xfactor, "order") <- order
    class(xfactor) <- c("smooth", class(xfactor))
    xfactor
   }
#-------------------------------------------------------------------------------------------
gamlss.ra <- function(x, y, w, df = NULL)
{
      #-----------------------------------------------------------------
      get.df <- function(lambda)
                    {  
                        aug <- sqrt(lambda)*aug
                       xaug <- as.matrix(rbind(dm,aug))
                         df <- sum(diag(solve(t(xaug) %*% (waug * xaug)) %*% (t(dm) %*% (waug[1:N] * dm))))
                         df
                     }
      #-----------------------------------------------------------------
       find.lambda.from.df <- function(df)
            {   
               usemode <- function(spar,df)
                    { 
                      ndf <- get.df(spar)
                      fun <- (ndf-df)**2
                      fun
                    }      
            a <- 0.000001 ;  c <- 10000; t <- 0.0001 ; r <- 0.61803399; z <- 1-r ; b <- r*a + z*c ;
            l4 <- a ; l3 <- c ; w1 <- 1 ;
            l1 <- if(c-b > b-a) b         else b-z*(b-a)
            l2 <- if(c-b > b-a) b+z*(c-b) else b 
            f1 <- usemode(l1,df)
            f2 <- usemode(l2,df)
            while(w1>0) 
                 { if(f2 < f1) { l4<-l1; l1<-l2 ; l2<- r*l1+z*l3; z8<-f1 ; f1<-f2 ; f2 <- usemode(l2,df) 
                               } 
                   else        { l3<-l2 ; l2<-l1; l1<- r*l2+z*l4; z8<-f2 ; f2<-f1 ; f1 <- usemode(l1,df) 
                               }
                   w1 <- abs(l3-l4) >  t*(abs(l1)+abs(l2))
                 }             
            lambda <- if(f1<f2) l1 else l2
            lambda
            }
       #---------------------------------------------------------------
       #---------------------------------------------------------------     
       find.lambda<-function(lambda)
           { 
           old.lambda <-lambda+1 # 1/(sum(coef(fit)^2+cov)/length(cov)) 
           niter <- 1
           cat("lambda=",lambda,"\n")
           while ( abs((old.lambda-lambda)/lambda) > 0.0001 && niter < 100 )
                { 
                 aug  <- sqrt(lambda)*aug
                xaug  <- as.matrix(rbind(dm,aug))
                 fit  <- lm.wfit(xaug-1,yaug,w=waug,method="qr") 
                 cov  <-  diag(chol2inv(fit$qr$qr[1:p, 1:p, drop = FALSE]))
           old.lambda <- lambda
               lambda <- 1/(sum(coef(fit)^2+cov)/length(cov))   
             cat("newlambda=",lambda,"\n")
                niter <- niter+1
                }
            lambda    
            }  
       #----------------------------------------------------------------
       # the ra() function starts here     
       dm <- as.matrix(attr(x,"design.matrix"))
      aug <- as.matrix(attr(x,"pen.augment"))
        N <- dim(dm)[1]
        p <- dim(dm)[2]
       aN <- dim(aug)[1]
       ap <- dim(aug)[2]
   lambda <- as.vector(attr(x,"lambda"))
       df <- as.vector(attr(x,"df"))
 estimate <- as.logical(attr(x,"estimate"))      
    if(p!=ap) stop("the dimensions of the augmented matrix and of the design matrix are incompatible")
    zeros <- rep(0,aN)
     ones <- rep(1,aN)
     yaug <- as.vector(c(y,zeros))
     waug <- as.vector(c(w,ones))
    # first if lambda not given find lambda from df 
    if(!is.null(df)&is.null(lambda)){ lambda <-  find.lambda.from.df(df)} # if df is set
    # second if the estimation of lambda is required
    if(estimate) { lambda<-find.lambda(lambda) } # if estimation of lambda
    # now given lambda fit the model
    if(!is.null(lambda))
        { 
     naug <- sqrt(lambda)*aug
     xaug <- as.matrix(rbind(dm,naug))
     #fit1 <- lsfit(xaug,yaug,wt=waug,intercept=F)
      if(is.null(attr(x,"explanatory")))
          { fit <- lm.wfit(xaug-1,yaug,w=waug,method="qr") 
            cov <-  diag(chol2inv(fit$qr$qr[1:p, 1:p, drop = FALSE]))
          nl.df <- sum(diag(solve(t(xaug) %*% (waug * xaug)) %*% (t(dm) %*% (waug[1:N] * dm))))
        coefSmo <- list(coef(fit), varcoeff=cov, lambda=lambda)
          }
     else
          {
           exdm <- as.matrix(attr(x,"explanatory"))
             eN <- dim(exdm)[1]
             ep <- dim(exdm)[2]
             if(eN!=aN) stop("the dimensions of the augmented matrix and of the design matrix for level 1 are incompatible") 
          Zeros <- matrix(rep(0,(N*ep)),ncol=ep)
          Zeros <- rbind(Zeros,exdm)
          bxaug <- cbind(xaug,Zeros) 
            fit <- lm.wfit(bxaug-1,yaug,w=waug,method="qr") 
            cov <-  diag(chol2inv(fit$qr$qr[1:p, 1:p, drop = FALSE])) #??
         # coef1 <- coef(fit)[(p+1):(p+ep)] # ??
          # cov1 <- diag(chol2inv(fit$qr$qr[(p+1):(p+ep), (p+1):(p+ep), drop = FALSE])) #??
          nl.df <- sum(diag(solve(t(xaug) %*% (waug * xaug)) %*% (t(dm) %*% (waug[1:N] * dm))))
          coefSmo <- list(coef(fit), varcoeff=cov,  lambda=lambda)
          }
      list(x = x, fitted.values = fitted(fit)[1:N], residuals = resid(fit)[1:N], 
         var = hat(sqrt(waug)*xaug,intercept=FALSE)[1:N], nl.df = nl.df, lambda = lambda, 
         coefSmo = coefSmo )
       }
}
