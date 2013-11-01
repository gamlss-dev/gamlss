# this is new and experimental to fit random coefficient models
# Mikis Stasinopoulos
rc <- function(formula, lambda=NULL) 
{
  rcParseFormula<-function (model) 
   {
     parseCond <- function(model) 
     {
        model <- eval(parse(text = paste("~", deparse(model))))[[2]]
        model.vars <- list()
        while (length(model) == 3 && (model[[1]] == as.name("*") || 
            model[[1]] == as.name("+"))) {
            model.vars <- c(model.vars, model[[3]])
            model <- model[[2]]
        }
        rev(c(model.vars, model))
     } # end of parseCond
     if (!inherits(model, "formula")) 
        stop("model must be a formula object")
     ans <- list(left = NULL, right = NULL, condition = NULL, left.name = character(0), 
            right.name = character(0))

 #   if (length(model) == 3) 
 #    {
 #       ans$left <- eval(model[[2]])
 #       if (inherits(ans$left, "POSIXt")) 
 #           ans$left <- as.POSIXct(ans$left)
 #       ans$left.name <- deparse(model[[2]])
 #    }
     model <- model[[length(model)]]
     if (length(model) == 3 && model[[1]] == as.name("|")) 
     {
        model.vars <- parseCond(model[[3]])
        ans$condition <- vector("list", length(model.vars))
        names(ans$condition) <- sapply(model.vars, deparse)
        for (i in seq(along = model.vars)) 
          {
            ans$condition[[i]] <- eval(model.vars[[i]],sys.frame(sys.parent(2))) #MSWednesday, April 16, 2003 at 09:58
            if (inherits(ans$condition[[i]], "POSIXt")) 
                ans$condition[[i]] <- as.POSIXct(ans$condition[[i]])
          }
        model <- model[[2]]
     }
     ans$right <- eval(model,sys.frame(sys.parent(2))) #MSWednesday, April 16, 2003 at 09:58
      if (inherits(ans$right, "POSIXt"))  ans$right <- as.POSIXct(ans$right)
       ans$right.name <- deparse(model)
    ans
   } # end of rcParseFormula
#---------------------------------------------------------------------- 
rcMatrix <- function(A,X)
      {
      for(ii in 1:dim(X)[1])
                { 
             P <- kronecker(A[ii,],X[ii,])
             if (ii==1) {PP <- P} 
             else {PP <- rbind(PP,P) }
      }
      PP
     }
#-----------------------------------------------------------------------
# For a design matrix of a grouping factor A create a matrix with zeros and ones
# indicating the row which the level of the factor is defined 
rcFactor <- function(A)
    {
      PP<-A[1,]
     
     for(ii in 2:dim(A)[1])
        {
         if ((any(A[(ii-1),]!=A[ii,]))) {P<-A[ii,]; PP <- rbind(PP,P)} 
        }
    PP
    }         
#-----------------------------------------------------------------------
# to be used by rc.getL
#   rc.getS<-function (object) 
#      {
#    Ncol <- round((-1 + sqrt(1 + 8 * length(object)))/2)
#    .C("matrixLog_pd", Factor = double(Ncol * Ncol), as.integer(Ncol), 
#        as.double(object), PACKAGE = "nlme")$Factor
#      }
#-----------------------------------------------------------------------
# from a positive definite matrix S this get L where L'L=S
# I am not sure about this any more? This is note used 
#   rc.getL <- function (object) 
#      {
#       require(nlme)
#        Ncol <- Dim(object)[2]
#        value <- array(rc.getS(object), c(Ncol, Ncol), attr(object, 
#            "Dimnames"))
#        attr(value, "logDet") <- sum(log(abs(svd(value)$d)))
#        value
#      }
#----------------------------------------------------------------------
# the function starts here 
    scall <- deparse(sys.call())
   # find out wht is in the formula     
    form <- rcParseFormula(formula)
   # the form$condition[[1]] should be the groups (factor)
    if(!inherits(form$condition[[1]], "factor")) stop("rc() expects a grouping factor as the conditional element of the formula")
   # get the Groups design matrix
   G<-model.matrix(~form$condition[[1]]-1,contrast="")
   # get the x-variables design matrix
   Xs <- model.matrix(~form$right)
   # the new design matrix with the x's in groups
   newX<-rcMatrix(G, Xs)
   # what it goes in the X matrix of the gamlss model  
    xfactor <- as.vector(form$right)  
   # checking the variancecovariance matrix
    if(is.null(lambda)) stop("lambda matrix not set")
    if(dim(lambda)[1]!=dim(lambda)[2]) stop("lambda matrix is not square") 
    if(dim(lambda)[1]!=dim(Xs)[2]) stop("the dimension of X not equal that of lambda")      
    aug <- solve(lambda)
    aug<-chol(aug)
    #aug <- rc.getL(aug)
    #aug <- eigen(aug)$vectors
    #now the bid aug matrix
    A<-rcFactor(G)
    Aug<-kronecker(A,aug)
    attr(xfactor, "design.matrix") <- Xs
    attr(xfactor, "big.design.matrix") <- newX 
    attr(xfactor, "groups") <- G  
    attr(xfactor, "pen.augment") <- aug
    attr(xfactor, "big.pen.augment") <- Aug
    attr(xfactor, "call") <- substitute(gamlss.rc(data[[scall]], z, w, df = df )) 
    attr(xfactor, "lambda") <- lambda
    class(xfactor) <- c("smooth", class(xfactor))
    xfactor
   }
#-------------------------------------------------------------------------------------------
 gamlss.rc <- function(x, y, w, df = NULL)
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
       # the gamlss.rc() function starts here     
       dm <- as.matrix(attr(x,"big.design.matrix"))
      aug <- as.matrix(attr(x,"big.pen.augment"))
        N <- dim(dm)[1]
        p <- dim(dm)[2]
       aN <- dim(aug)[1]
       ap <- dim(aug)[2]
   lambda <- as.matrix(attr(x,"lambda"))
    if(p!=ap) stop("the dimensions of the augmented matrix and of the design matrix are incompatible")
    zeros <- rep(0,aN)
     ones <- rep(1,aN)
     yaug <- as.vector(c(y,zeros))
     waug <- as.vector(c(w,ones))
     xaug <- as.matrix(rbind(dm,aug))
     #fit1 <- lsfit(xaug,yaug,wt=waug,intercept=F)
      fit <- lm.wfit(xaug-1,yaug,w=waug,method="qr") 
      cov <-  diag(chol2inv(fit$qr$qr[1:p, 1:p, drop = FALSE]))
    nl.df <- sum(diag(solve(t(xaug) %*% (waug * xaug)) %*% (t(dm) %*% (waug[1:N] * dm))))
    # equivalent way to estimate df Hodges ansd Sargent (1998)
    #   s<-svd(sqrt(waug)*xaug)
    #   U1<-s$u[1:N,]
    #  df<- sum(diag(U1 %*% t(U1) ) )
    
    #test<- lm(fitted(fit)[1:N]~x,weights=w) 
    
     list(x = x, fitted.values=fitted(fit)[1:N], residuals=resid(fit)[1:N], #residuals=y-resid(test)
         var=hat(sqrt(waug)*xaug,intercept=FALSE)[1:N], nl.df =nl.df, lambda=lambda, 
         coefSmo=list(coef(fit), varcoeff=cov) )
    }
#-------------------------------------------------------------------------------------------
