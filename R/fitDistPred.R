# Created by Mikis Stasinopoulos February 2018
# based on the gamlssVGD() function
# this function fits a gamlssML model but also save the prediction deviance
# it creates  gamlssML and gamlssVGD objects
# TO DO
#  i) the prediction residuals are not working for discrete data (probably neither in gamlssVGD())
#  ii) 
#-----------------------------------------------------------------------------
gamlssMLpred <-function(response = NULL, 
                            data = NULL, # original data 
                          family = NO,  
                            rand = NULL, # 1 for training 2 for validation
                         newdata = NULL, 
                                 ...)
{
#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
# local function -------------------------------------------------------------  
predictAll.V <- function(object, newdata=NULL, se.fit=FALSE)
  {
     y  <- if (is.null(newdata)) y else newdata
    out <- list(y=y)
    if ("mu" %in% object$par)   
      out$mu <- predict(object, newdata=newdata, what = "mu", se.fit = se.fit )
    if ("sigma" %in% object$par)  
      out$sigma <- predict(object,  newdata=newdata, what = "sigma", se.fit = se.fit)
    if ("nu" %in% object$par)     
      out$nu <- predict(object,  newdata=newdata, what = "nu",  se.fit = se.fit)
    if ( "tau" %in% object$par)    
      out$tau <- predict(object,  newdata=newdata, what = "tau", se.fit = se.fit)
    attr(out, "family") <- object$family
    return(out)
}
#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
# this is to replicate rqres within gamlssMLpred enviroment
# it is used as in gamlss()
rqres <- function (pfun = "pNO", 
                   type = c("Continuous", "Discrete", "Mixed"),
                   censored = NULL,  
                   ymin = NULL, 
                   mass.p = NULL, 
                   prob.mp = NULL,
                   y = y,
                   ... )
{ }
body(rqres) <-  eval(quote(body(rqres)), envir = getNamespace("gamlss"))
##---------------------------------------------------------------------------------------
##---------------------------------------------------------------------------------------  
#-------------------------------------------------------------------------------
# main function starts here
#-------------------------------------------------------------------------------
  if (is.null(rand)&&is.null(newdata)) stop("rand or newdata should be set")
  if (!is.null(rand)&&!is.null(newdata)) stop("only rand or newdata should be set NOT both")
  if (!is.null(rand))
  {
    if ( any(!rand%in%c(1,2))) stop("rand values should be 1 or 2")
    resTraining <- subset(response, rand==1)
       resValid <- subset(response, rand==2)   
  }  
  fname <- as.gamlss.family(family)
   dfun <- paste("d", fname$family[[1]],sep="")
   pfun <- paste("p", fname$family[[1]],sep="")
   lpar <- length(fname$parameters)
  dtype <- fname$type
  # FIT MODEL + predict --------------------------------------------------------
  if (!is.null(rand))#  if rand is set to this ---------------------------------
  {
          m1 <- gamlssML(resTraining, family=fname,  ...) 
     nfitted <- predictAll.V(m1, newdata=resValid)
dim1newdata  <- length(resValid)
  } else
  {
          m1 <- gamlssML(response, family=fname,  ...) 
     nfitted <- predictAll.V(m1, newdata=newdata)
dim1newdata  <- length(newdata)   
  }  
#------------------------------------------------------------------------------
# get  y for new data  if binomial       
if (fname$family[1] %in% .gamlss.bi.list)# if binomial
    {
      if (NCOL(nfitted$y) == 1) 
      {
        y1 <- nfitted$y 
      }
      else                 
      {
        bd <- nfitted$y[,1] + nfitted$y[,2]
        y1 <- nfitted$y[,1]
      }
    } else
    {
      y1 <- nfitted$y 
    }
#---------------------------------------------------------------------
# jump depending on the number of parameters 
    if(lpar==1) 
    {
      if (fname$family[[1]] %in% .gamlss.bi.list)
      {
        devi <- call(dfun, x=y1, mu =  nfitted$mu, bd=bd, log=TRUE) 
        ures <-  call("rqres", pfun=pfun, type=dtype,
                      ymin=fname$rqres[[1]][["ymin"]], 
                      y=y1, bd=bd, mu= nfitted$mu)   
      } else
      {
        devi <- call(dfun, x=y1, mu =  nfitted$mu,  log=TRUE) 
        ures <- call("rqres", pfun=pfun, type=dtype, 
                      ymin=fname$rqres[[1]][["ymin"]], 
                      y=y1, mu= nfitted$mu)   
      } 
    }
    else if(lpar==2)
    {
      if (fname$family[[1]] %in% .gamlss.bi.list)
      {
        devi <- call(dfun, x=y1, mu =  nfitted$mu, sigma= nfitted$sigma,  
                     bd=bd, log=TRUE) 
        ures <-  call("rqres", pfun=pfun, type=dtype, 
                      ymin=fname$rqres[[1]][["ymin"]], 
                      y=y1, mu= nfitted$mu, sigma= nfitted$sigma, bd=bd)
      } else
      {
        devi <- call(dfun, x=y1, mu =  nfitted$mu, sigma =  nfitted$sigma, log=TRUE) 
        ures <- call("rqres", pfun=pfun,  type=dtype, 
                     ymin=fname$rqres[[1]][["ymin"]], 
                     y=y1, mu= nfitted$mu, sigma= nfitted$sigma )
      } 
    }
    else if(lpar==3)
    {
      if (fname$family[[1]] %in% .gamlss.bi.list)
      {
        devi <-  call(dfun, x=y1, mu =  nfitted$mu, sigma =  nfitted$sigma, 
                      nu =  nfitted$nu, bd=bd, log=TRUE)
        ures <- call("rqres", pfun=pfun,  type=dtype, 
                     ymin=fname$rqres[[1]][["ymin"]], 
                     y=y1, mu= nfitted$mu, sigma= nfitted$sigma,
                     nu =  nfitted$nu, bd=bd)
      } else
      {
        devi <-  call(dfun, x=y1, mu = nfitted$mu, sigma =  nfitted$sigma, 
                      nu =  nfitted$nu, log=TRUE)
        ures <- call("rqres", pfun=pfun,  type=dtype, 
                     ymin=fname$rqres[[1]][["ymin"]], 
                     y=y1, mu= nfitted$mu, sigma= nfitted$sigma,
                     nu =  nfitted$nu)
      } 
    }
    else 
    {
      if (fname$family[[1]] %in% .gamlss.bi.list)
      {
        devi <-  call(dfun, x=y1, mu = nfitted$mu, sigma =  nfitted$sigma, 
                      nu =  nfitted$nu, tau = nfitted$tau, bd=bd, log=TRUE)
        ures <- call("rqres", pfun=pfun,  type=dtype, 
                     ymin=fname$rqres[[1]][["ymin"]], 
                     y=y1, mu= nfitted$mu, sigma= nfitted$sigma,
                     nu =  nfitted$nu, tau = nfitted$tau, bd=bd)
      } else
      {
        devi <-  call(dfun, x=y1, mu =  nfitted$mu, sigma =  nfitted$sigma,
                      nu =  nfitted$nu, tau = nfitted$tau, log=TRUE)
        ures <- call("rqres", pfun=pfun,  type=dtype, 
                     ymin=fname$rqres[[1]][["ymin"]], 
                     y=y1, mu= nfitted$mu, sigma= nfitted$sigma,
                     nu =  nfitted$nu, tau= nfitted$tau)
      } 
    }
         Vresid <- eval(ures)
       dev.incr <- -2 * eval(devi)
            dev <- sum(dev.incr)
         m1$VGD <- dev
     m1$IncrVGD <- dev.incr
m1$predictError <- dev/dim1newdata
    m1$residVal <- Vresid 
  #   }#----------END of newdata -------------------------------------------------
  class(m1) <- c("gamlssVGD", "gamlss", "gam",    "glm",    "lm"  )
  m1
} # end of function
#-----------------------------------------------------------------------------
#-------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------
fitDistPred <- function(y,
                    type = c("realAll", "realline", "realplus","real0to1","counts", "binom" ), 
              try.gamlss = FALSE,  # whether to try the gamlss() if gamlssML() fails
                   extra = NULL,  # for extra distributions to include 
                    data = NULL, 
                    rand = NULL, # 1 for training 2 for validation
                 newdata = NULL,
                   trace = FALSE,
                    ...)
{
  # if (!is.null(data)) {attach(data); on.exit(detach(data))}
  #if (!is.null(data)) {attach(data, name="TheDatA"); on.exit(detach(TheDatA))}
  if (is.null(rand)&&is.null(newdata)) stop("rand or newdata should be set")
  if (!is.null(rand)&&!is.null(newdata)) stop("only rand or newdata should be set NOT both")
     y <- if (!is.null(data)) get(deparse(substitute(y)), envir=as.environment(data)) else y
  type <- match.arg(type)
  DIST <- switch(type, "realAll"=.realAll, 
                      "realline"=.realline, 
                      "realplus"=.realplus,
                      "real0to1"=.real0to1,
                        "counts"=.counts,
                         "binom"=.binom 
  )
  if  (!is.null(extra)) DIST <- unique(c(DIST, extra))
  # do we need weights here 
  m0 <- switch(type,  "realAll"= gamlssMLpred(y, rand=rand, newdata=newdata, family=NO, ...),
                     "realline"= gamlssMLpred(y, rand=rand, newdata=newdata, family=NO, ...), 
                     "realplus"= gamlssMLpred(y, rand=rand, newdata=newdata,family=EXP, ...),
                     "real0to1"= gamlssMLpred(y, rand=rand, newdata=newdata,family=BE, ...),
                       "counts"= gamlssMLpred(y, rand=rand, newdata=newdata,family=PO, ...),
                        "binom"= gamlssMLpred(y, rand=rand, newdata=newdata,family=BI, ...) 
  ) 
  failed <- list() 
    fits <- list()
      pb <- txtProgressBar(max = length(DIST), style=3)
  for (i in 1:length(DIST)) 
  {
    setTxtProgressBar(pb, i)    
    m1 <- try(gamlssMLpred(y, rand=rand, newdata=newdata,family=DIST[i], ...), silent=TRUE)
    if (any(class(m1)%in%"try-error")&&try.gamlss==TRUE) 
    { 
      m1 <-  try(gamlssVGD(y~1,family=DIST[i], trace=FALSE, ...),  silent=TRUE)
    }
    
    if (any(class(m1)%in%"try-error"))
    {
      failed <- c(failed, DIST[i]) 
    }
    else
    {
      vgd <- VGD(m1)
      names(vgd) <- DIST[i]
 if (trace)     cat(DIST[i], vgd, "\n")
      fits <- c(fits, vgd)
      if ( VGD(m1) <  VGD(m0)) 
      {
        m0<-m1 
      }
    }
  }
      close(pb)        
  m0$failed <- failed
  fits <- unlist(fits)
  m0$fits <- fits[order(fits)]         
  m0  
}
#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
