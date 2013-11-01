# Authors Mikis Stasinopoulos Bob Rigby and Vlasios Voudouris
# created 11-04-12
#--------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------
#This function is design to help the user to construct centile estimation.
#It is only applicable when only "one" explanatory variable is available (usually age).
#It stats by fitting a normal error distribution and smooth function for mu and then proceeds by fitting #several appropriate distributions.
#The set of gamlss.family distribution to fit are specified in the argument families. 
#The default families arguments is LMS=c("BCCGo",  "BCPEo", "BCTo") that is the LMS class of distributions.
#Note that this class is only appropriate when y is positive (with no zeros). If the response variable contains negative values and/or zeros then use 
#the argument theSHASH theSHASH <-  c("NO", "SHASHo") or add any other distribution which you think is appropriate
LMS <- c("BCCGo",  "BCPEo", "BCTo")
theSHASH <-  c("NO", "SHASHo")
#------------
lms <- function(y, x,
        families = LMS,
            data = NULL, 
               k = 2, # for the AIC
            cent = 100*pnorm((-4:4)*2/3),
     calibration = TRUE,
         trans.x = FALSE,  
       lim.trans = c(0.001, 1.5),        
          legend = FALSE,
           mu.df = NULL,
        sigma.df = NULL,
           nu.df = NULL,
          tau.df = NULL,
       method.pb = c("ML", "GAIC"),
              ... 
                )  
{
##  require(gamlss)
## the families to fit
        FAM <- families      
## which method
  method.pb <- match.arg(method.pb)
## get the variables  
       ylab <- deparse(substitute(y))
       xlab <- deparse(substitute(x))
          y <- if (!is.null(data)) get(deparse(substitute(y)), envir=as.environment(data)) else y
          x <- if (!is.null(data)) get(deparse(substitute(x)), envir=as.environment(data)) else x
## ----------------------  
## starting value model (we assume that this will work). Note no sigma is fitted here
## if need to chech for transformation in x        
  if (trans.x) # if x^p
  {
    cat("*** Checking for transformation for x ***", "\n") 
    switch(method.pb, 
      "ML"= {  fn <- function(p) deviance(gamlss(y~pb(x^p), data=data, c.crit = 0.01, trace=FALSE))
              par <- optim(.5, fn, lower=lim.trans[1], upper=lim.trans[2], method="L-BFGS-B")$par
               nx <- x^par 
               ox <- x
                x <- x^par 
        #       m0 <- gamlss(y~pb(x), data=data, c.crit = 0.01)
            },
    "GAIC"= {  fn <- function(p) deviance(gamlss(y~pb(x^p, method="GAIC", k=k), data=data, c.crit = 0.01, trace=FALSE))
              par <- optim(.5, fn, lower=lim.trans[1], upper=lim.trans[2], method="L-BFGS-B")$par
               ox <- x
                x <- x^par 
        #      m0 <- gamlss(y~pb(x, method="GAIC", k=k ), data=data, c.crit = 0.01)
            })   
    cat('*** power parameters ', par,"***"," \n") 
  } 
## now fit the model  
    cat('*** Initial  fit***'," \n")       
    switch(method.pb, 
        "ML"= {m0 <- gamlss(y~pb(x), sigma.formula=~1, data=data, c.crit = 0.01)},
      "GAIC"= {m0 <- gamlss(y~pb(x, method="GAIC", k=k), sigma.formula=~1, data=data, c.crit = 0.01)}) ## initial fit  finish
## creating lists etc 
     failed <- list() 
       fits <- list()
        aic <- AIC(m0, k=k)
       fits <- c(fits, aic) 
## fitting the diferent models in FAM        
  for (i in 1:length(FAM)) 
  {
    cat('*** Fitting', FAM[i], "***","\n")  
     switch(method.pb, 
         "ML"= { m1 <- try(gamlss(y ~ pb(x, df=mu.df),
              sigma.formula = ~pb(x, df=sigma.df),
                 nu.formula = ~pb(x, df=nu.df), 
                tau.formula = ~pb(x, df=tau.df), 
                family = FAM[i], data = data,
               mu.start = fitted(m0), ...), silent=TRUE)},
        "GAIC"= { m1 <- try(gamlss(y ~ pb(x,  method="GAIC", k=k, df=mu.df),
              sigma.formula = ~pb(x,  method="GAIC", k=k, df=sigma.df),
                 nu.formula = ~pb(x,  method="GAIC", k=k, df=nu.df), 
                tau.formula = ~pb(x,  method="GAIC", k=k, df=tau.df), 
                family = FAM[i], data = data,
               mu.start = fitted(m0), ...), silent=TRUE)
         })      
    if (any(class(m1)%in%"try-error")) # if fitting failed
    {
      cat(FAM[i], " failed", "\n")
          failed <- c(failed, FAM[i]) 
    }
    else
    {
             aic <- AIC(m1, k=k)
      names(aic) <- FAM[i]
            fits <- c(fits, aic)
      if (AIC(m1, k=k) < AIC(m0, k=k)) 
      {
        m0<-m1 
      }
    }
  }
## transformation needed         
     if (trans.x)   
       { 
          x <- ox
   m0$power <- par 
       }
## save the rest information        
  m0$failed <- failed
       fits <- unlist(fits)
    m0$fits <- fits[order(fits)] 
    m0$xvar <- x#with(DaTa,x)
       m0$y <- y#with(DaTa,y)
  if (!is.null(data)) m0$call$data  <- substitute(data)
## calibration
  if (calibration)
  {
    calibration(m0, xvar=x, cent=cent, pch = 15, cex = 0.5, col = gray(0.7), ylab=ylab, xlab=xlab, legend=legend)	
  } 
  else 
  {
    centiles(m0, xvar=x, cent=cent, pch = 15, cex = 0.5, 
             col = gray(0.7), ylab=ylab, xlab=xlab, legend=legend)		
  }
  m0  # save the last model
}
#---------------------------------------------------------------------------------
#---------------------------------------------------------------------------------
#---------------------------------------------------------------------------------
# this function is appropriate to used when fitted model fails to c
calibration <- function(object, xvar, cent=100*pnorm((-4:4)*2/3), legend=FALSE, fan=FALSE,  ...)
{
  z   <-  quantile(resid(object), probs = cent/100)
  p   <-  pNO(z, mu=0, sigma=1)
  percent <- 100*p
  if (fan)
  {
    centiles.fan(object, xvar=xvar, cent=percent,   ...)  
  }
  else
  {
    centiles(object, xvar=xvar, cent=percent, legend=legend,  ...)
  }
}
#---------------------------------------------------------------------------------
#---------------------------------------------------------------------------------

#---------------------------------------------------------------------------------
#---------------------------------------------------------------------------------
#test <- function(y,x, data)
#{
#  #Data<-assign(deparse(substitute(data)), get(deparse(substitute(data))))
#  #currentEnv <- environment()   
#  y <- if (!is.null(data)) get(deparse(substitute(y)), envir=as.environment(data)) else y
#  x <- if (!is.null(data)) get(deparse(substitute(x)), envir=as.environment(data)) else x
# # assign(deparse(substitute(x)),get(deparse(substitute(x)), envir=as.environment(data)))
##  assign(deparse(substitute(y)),get(deparse(substitute(y)), envir=as.environment(data)))
# # browser()
#  m1 <- gamlss(y~x)
# # with(Data, {m1<- gamlss(y~x); assign("m1", m1, envir=currentEnv)})
#  m1
#}
# get the data
#if (!is.null(data)) 
#  {
#   assign(deparse(substitute(x)),get(deparse(substitute(x)), envir=as.environment(data)))
#   assign(deparse(substitute(y)),get(deparse(substitute(y)), envir=as.environment(data)))
#          DaTa <- data.frame(y=get(deparse(substitute(y)), envir=as.environment(data)), x=get(deparse##(substitute(x)), envir=as.environment(data)))3
#  names(DaTa) <- c(deparse(substitute(y)), deparse(substitute(x)))
#  DaTa<-assign(deparse(substitute(data)), get(deparse(substitute(data))))
#  }
#if (is.null(data))  
#  {}
#--------------------------------------  