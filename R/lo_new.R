##---------------------------------------------------------------------------------------
## the lo(), lo.control() and gamlss.lo() functions
## are based on R loess() and S-plus lo() 
## author Mikis Stasinopoulos
##---------------------------------------------------------------------------------------
## created by MS : 13-08-12
## based on Brian Ripley impementation of loess in R function
## This replace the previous version 2002 version
#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
lo <-function(formula, control=lo.control(...), ...) 
{ 
#------------------------------------------
# function starts here
#------------------------------------------
    scall <- deparse(sys.call())
if (!is(formula, "formula")) stop("lo() needs a formula starting with ~")
# get where "gamlss" is in system call
# it can be in gamlss() or predict.gamlss()  
    rexpr <- grepl("gamlss",sys.calls()) ## 
for (i in length(rexpr):1)
   { 
 position <- i # get the position
 if (rexpr[i]==TRUE) break
   }
  # 
gamlss.env <- sys.frame(position) #gamlss or predict.gamlss
##---
## get the data
if (sys.call(position)[1]=="predict.gamlss()")
     { # if predict is used 
      Data <- get("data", envir=gamlss.env)
     }
else { # if gamlss() is used
	#stop("the option data in gamlss() is required for lo() to work")
     if (is.null(get("gamlsscall", envir=gamlss.env)$data)) 
         { # if no data argument but the formula can be interpreted
     Data <- model.frame(formula)	
         }
     else
         {# data argument in gamlss 
     Data <- get("gamlsscall", envir=gamlss.env)$data
         }
     }
     Data <- data.frame(eval(substitute(Data)))
     #===== 
      len <- dim(Data)[1] # get the lenth of the data
    #  get the free arguments 
    alist <- list(...)
    # check if df are defined
    if (!is.null(alist$df))  control$enp.target <- alist$df 
## out
     xvar <- rep(0,  len) # model.matrix(as.formula(paste(paste("~",ff, sep=""), "-1", sep="")), data=Data) #
   attr(xvar,"formula")     <- formula
    attr(xvar,"control")    <- control 
   attr(xvar, "gamlss.env") <- gamlss.env
   attr(xvar, "data")       <- as.data.frame(Data)
   attr(xvar, "call")       <- substitute(gamlss.lo(data[[scall]], z, w, ...)) 
   attr(xvar, "class")      <- "smooth"
   xvar
}
#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
lo.control <-  function (span = 0.75, enp.target=NULL, degree = 2,
      parametric = FALSE, drop.square = FALSE, normalize = TRUE,
      family = c("gaussian", "symmetric"),
      method = c("loess", "model.frame"),
      surface = c("interpolate", "direct"), 
      statistics = c("approximate", "exact"), 
      trace.hat = c("exact", "approximate"), 
      cell = 0.2, 
      iterations = 4, ...) 
{
       list(span = span, enp.target=enp.target, degree=degree,
         parametric = parametric, drop.square = drop.square, normalize = normalize, family=family,
         method=method, surface = surface, statistics = statistics, trace.hat = trace.hat, 
         cell = cell, iterations = iterations)
}
#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
# the definition of the backfitting additive function
gamlss.lo <-function(x, y, w, xeval = NULL, ...)
{           
    formula <- attr(x,"formula")
    formula <- as.formula(paste("Y.var",deparse(formula), sep=""))
    control <- as.list(attr(x, "control"))
#gamlss.env <- as.environment(attr(x, "gamlss.env"))
      OData <- attr(x,"data") 
       Data <-  if (is.null(xeval)) OData #the trick is for prediction
               else  OData[seq(1,length(y)),]
       Y.var <- y
       W.var <- w     
      Data <- data.frame(eval(substitute(Data)),Y.var,W.var)   
       if (is.null(control$enp.target))
        { 
       fit <- loess(formula, data=Data, weights=W.var, 
                     span=control$span, 
                     degree=control$degree, normalize=control$normalize,
                     family=control$family, 
                     control=loess.control(surface=control$surface,
                     statistics=control$statistics,
                     trace.hat=control$trace.hat,
                      iterations= control$iterations)) 
         }
       else 
         { 
       	fit <-loess(formula, data=Data, weights=W.var,
                     enp.target=control$enp.target, 
                     degree=control$degree, normalize=control$normalize,
                     family=control$family, 
                     control=loess.control(surface=control$surface,
                     statistics=control$statistics,
                     trace.hat=control$trace.hat,
                      iterations= control$iterations)) 
         } 
        df <- fit$trace.hat-1 
        fv <- fitted(fit) 
 residuals <- Y.var-fv
  if (is.null(xeval))
    {
   list(fitted.values=fv, residuals=residuals,
     nl.df = df, lambda=fit$pars[["span"]], ## we nead df's here 
     coefSmo = fit, var=NA)    # var=fv has to fixed
    }
else 
    {
   ll<-dim(OData)[1]
   pred <- predict(fit,newdata = OData[seq(length(y)+1,ll),])
    }         
}
#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------      
