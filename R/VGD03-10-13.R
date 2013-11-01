#----------------------------------------------------------------------------------------
#rand <- sample(2, 610, replace=T, prob=c(0.6,0.4))
#----------------------------------------------------------------------------------------
# this is the original function
VGD <-function(formula = NULL, 
        sigma.formula =~1, 
           nu.formula =~1, 
          tau.formula =~1, 
                 data = NULL,
               family = NO,  
              control = gamlss.control(trace=FALSE),
                 rand = NULL,
                 ...)
 {
       fname <- as.gamlss.family(family)
        dfun <- paste("d", fname$family[[1]],sep="")
        lpar <- length(fname$parameters)
        if (is.null(formula)) stop("no formula is set in VGD")
        if (is.null(rand)) stop("no subset is set in VGD")
        if ( any(!rand%in%c(1,2))) stop("rand values should ne 1 or 2")
        if (is.null(data)) stop("the data argument is needed in VGD")
      dataor <- subset(data, rand==1)
      datava <- subset(data, rand==2)    
          m1 <- gamlss(formula=formula, sigma.formula=sigma.formula, nu.formula=nu.formula, tau.formula=tau.formula,
                      data=dataor, family=family, control=control, ...) #
         nmu <- predict(m1,newdata=datava, type="response", data=dataor)
     if ("sigma"%in% names(fname$par))  
      nsigma <- predict(m1,newdata=datava, type="response", data=dataor,what="sigma")
     if (  "nu" %in% names(fname$par))     
         nnu <- predict(m1,newdata=datava, type="response", data=dataor,what="nu")
     if ( "tau" %in% names(fname$par))    
        ntau <- predict(m1,newdata=datava, type="response", data=dataor,what="tau")
         ny <-  model.extract(model.frame(formula,data=datava),"response")
        # ny <- datava$y # corrected Monday, November 6, 2006 DS
      if (fname$family[[1]] %in% .gamlss.bi.list)# if binomial
              {
               if (NCOL(ny) == 1) 
                 {
                   y1 <- ny 
                 }
                else                 
                 {
                 bd <- ny[, 1] + ny[, 2]
                 y1 <- ny[, 1]
                  }
              }
      if(lpar==1) 
       {
       devi <- if (fname$family[[1]] %in% .gamlss.bi.list) call(dfun, x=y1, mu = nmu, bd=bd, log=TRUE) 
                else call(dfun, x= ny, mu = nmu, log=TRUE) 
       }
      else if(lpar==2)
       {
       devi <- if (fname$family[[1]] %in% .gamlss.bi.list) call(dfun, x=y1, mu = nmu, sigma=nsigma,  bd=bd, log=TRUE) 
                else call(dfun, x=ny, mu = nmu, sigma = nsigma, log=TRUE) 
       }
      else if(lpar==3)
       {
       devi <-  if (fname$family[[1]] %in% .gamlss.bi.list) call(dfun, x=y1, mu = nmu, sigma = nsigma, nu = nnu, bd=bd, log=TRUE)
                else  call(dfun, x=ny, mu = nmu, sigma = nsigma, nu = nnu, log=TRUE)
       }
      else 
       {
       devi <-  if (fname$family[[1]] %in% .gamlss.bi.list) call(dfun, x=y1, mu = nmu, sigma = nsigma, nu = nnu, tau=ntau, bd=bd, log=TRUE)
                else  call(dfun, x=ny, mu = nmu, sigma = nsigma,nu = nnu,tau = ntau, log=TRUE)
       
       }
       dev <- -2*sum(eval(devi))
     dev
 }
#----------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------
#
#----------------------------------------------------------------------------------------
# identical to VGD put the output is defferent 
VGD1 <-function(formula = NULL, 
        sigma.formula =~1, 
           nu.formula =~1, 
          tau.formula =~1, 
                 data = NULL,
               family = NO,  
              control = gamlss.control(trace=FALSE),
                 rand = NULL,
                 ...)
 {
       fname <- as.gamlss.family(family)
        dfun <- paste("d", fname$family[[1]],sep="")
        lpar <- length(fname$parameters)
        if (is.null(formula)) stop("no formula is set in VGD")
        if (is.null(rand)) stop("no subset is set in VGD")
        if ( any(!rand%in%c(1,2))) stop("rand values should ne 1 or 2")
        if (is.null(data)) stop("the data argument is needed in VGD")
      dataor <- subset(data, rand==1)
      datava <- subset(data, rand==2)    
          m1 <- gamlss(formula=formula, sigma.formula=sigma.formula, nu.formula=nu.formula, tau.formula=tau.formula,
                      data=dataor, family=family, control=control, ...) #
         fit <- predictAll(m1,newdata=datava, type="response", data=dataor)
         if (fname$family[[1]] %in% .gamlss.bi.list)# if binomial
              {
               if (NCOL(fit$y) == 1) 
                 {
                   y1 <- fit$y 
                 }
                else                 
                 {
                 bd <- fit$y[,1] + fit$y[,2]
                 y1 <- fit$y[,1]
                 }
              }
      if(lpar==1) 
       {
        devi <- if (fname$family[[1]] %in% .gamlss.bi.list) call(dfun, x=y1, mu = fit$mu, bd=bd, log=TRUE) 
                else call(dfun, x= fit$y, mu = fit$mu, log=TRUE)  
       }
      else if(lpar==2)
       {
       devi <- if (fname$family[[1]] %in% .gamlss.bi.list) call(dfun, x=y1, mu = fit$mu, sigma=fit$sigma, bd=bd, log=TRUE) 
                else call(dfun, x=fit$y, mu = fit$mu, sigma = fit$sigma, log=TRUE) 
       }
      else if(lpar==3)
       {
       devi <- if (fname$family[[1]] %in% .gamlss.bi.list) call(dfun, x=y1, mu = fit$mu, sigma=fit$sigma, nu = fit$nu, bd=bd, log=TRUE) 
                else call(dfun, x=fit$y, mu = fit$mu, sigma = fit$sigma, nu = fit$nu, log=TRUE)  
       }
      else 
       {
       devi <-  if (fname$family[[1]] %in% .gamlss.bi.list) call(dfun, x=y1, mu = fit$mu, sigma=fit$sigma, nu = fit$nu, tau = fit$tau, bd=bd, log=TRUE) 
                else call(dfun, x=fit$y, mu = fit$mu, sigma = fit$sigma,nu = fit$nu,tau = fit$tau, log=TRUE)
       }
       dev <- -2*sum(eval(devi))
     ll<-list(oldGD=deviance(m1), newGD=dev, oldPE=deviance(m1)/m1$noObs, newPE=dev/dim(datava)[1])
     ll
 }
#----------------------------------------------------------------------------------------
# identical to VGD1 but uses newdata rather that rand as argument
VGD2<- function (formula = NULL, 
           sigma.formula = ~1, 
              nu.formula = ~1, 
             tau.formula = ~1, 
                    data = NULL, 
                  family = NO, 
                 control = gamlss.control(trace = FALSE), 
                 newdata = NULL, 
                 ...) 
{
    fname <- as.gamlss.family(family)
    dfun <- paste("d", fname$family[[1]], sep = "")
    lpar <- length(fname$parameters)
    if (is.null(formula)) 
        stop("no formula is set in VGD")
    if (is.null(newdata)) 
        stop("no neadata is set in VGD")
    if (is.null(data)) 
        stop("the data argument is needed in VGD")
    m1 <- gamlss(formula = formula, sigma.formula = sigma.formula, nu.formula = nu.formula, 
        tau.formula = tau.formula, data = data, family = family, control = control, 
        ...)
    nfitted <- predictAll(m1, newdata=newdata, data=data)
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
              }
    
    if (lpar == 1) 
    {
    devi <- if (fname$family[[1]] %in% .gamlss.bi.list) call(dfun, x=y1, mu = nfitted$mu, bd=bd, log=TRUE) 
                else call(dfun, x = nfitted$y, mu = nfitted$mu, log = TRUE)
    }
    else if (lpar == 2) 
    {
    devi <- if (fname$family[[1]] %in% .gamlss.bi.list) call(dfun, x=y1, mu = nfitted$mu, sigma=nfitted$sigma, bd=bd, log=TRUE) 
                else call(dfun, x = nfitted$y, mu =  nfitted$mu, sigma = nfitted$sigma, log = TRUE)
    }
    else if (lpar == 3) 
    {
    devi <- if (fname$family[[1]] %in% .gamlss.bi.list) call(dfun, x=y1, mu = nfitted$mu, sigma=nfitted$sigma, nu = nfitted$nu, bd=bd, log=TRUE) 
                else call(dfun, x = nfitted$y, mu =  nfitted$mu, sigma = nfitted$sigma, nu = nfitted$nu, log = TRUE)
    }
    else 
    {
    devi <- if (fname$family[[1]] %in% .gamlss.bi.list) call(dfun, x=y1, mu = nfitted$mu, sigma=nfitted$sigma, nu = nfitted$nu, tau = nfitted$tau, bd=bd, log=TRUE) 
                else call(dfun, x = nfitted$y, mu =  nfitted$mu, sigma = nfitted$sigma, nu = nfitted$nu, tau = nfitted$tau, log = TRUE)
    }
    dev <- -2 * sum(eval(devi))
    ll<-list(oldGD=deviance(m1), newGD=dev, oldPE=deviance(m1)/m1$noObs, newPE=dev/dim(newdata)[1])
    ll
}
#----------------------------------------------------------------------------------------
# this requires a fitted model and the new data only
TGD<- function (object,   newdata = NULL, ...) 
{
   if (!is.gamlss(object)) stop("not a gamlss object")
    fname <- as.gamlss.family(object$family[[1]])
    dfun <- paste("d", fname$family[[1]], sep = "")
    lpar <- length(fname$parameters)
    if (is.null(newdata)) 
        stop("no newdata is set in VGD")
   # if (is.null(data)) 
   #     stop("the data argument is needed in VGD")
    nfitted <- predictAll(object, newdata=newdata, ...)
    if (is.null(nfitted$y)) stop("the response variables is missing in the newdata")
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
              }
    if (lpar == 1) 
    {
    devi <- if (fname$family[[1]] %in% .gamlss.bi.list) call(dfun, x=y1, mu = nfitted$mu, bd=bd, log=TRUE) 
                else call(dfun, x = nfitted$y, mu = nfitted$mu, log = TRUE)
                     call(dfun, x = nfitted$y, mu = nfitted$mu, log = TRUE)
    }
    else if (lpar == 2) 
    {
    devi <- if (fname$family[[1]] %in% .gamlss.bi.list) call(dfun, x=y1, mu = nfitted$mu, sigma=nfitted$sigma, bd=bd, log=TRUE) 
                else call(dfun, x = nfitted$y, mu =  nfitted$mu, sigma = nfitted$sigma, log = TRUE)
    }
    else if (lpar == 3) 
    {
    devi <- if (fname$family[[1]] %in% .gamlss.bi.list) call(dfun, x=y1, mu = nfitted$mu, sigma=nfitted$sigma,  nu = nfitted$nu, bd=bd, log=TRUE) 
                else  call(dfun, x = nfitted$y, mu =  nfitted$mu, sigma = nfitted$sigma, nu = nfitted$nu, log = TRUE)
    }
    else 
    {
    devi <- if (fname$family[[1]] %in% .gamlss.bi.list) call(dfun, x=y1, mu = nfitted$mu, sigma=nfitted$sigma,  nu = nfitted$nu, tau = nfitted$tau, bd=bd, log=TRUE) 
                else  call(dfun, x = nfitted$y, mu =  nfitted$mu, sigma = nfitted$sigma, nu = nfitted$nu, tau = nfitted$tau, log = TRUE)
    }
    dev <- -2 * sum(eval(devi))
    ll<-list(oldGD=deviance(object), newGD=dev, oldPE=deviance(object)/object$noObs, newPE=dev/dim(newdata)[1])
    ll
}
#----------------------------------------------------------------------------------------    
     
