################################################################################
################################################################################
################################################################################
################################################################################
# this work is a rethinking on Training, Validation and Test data sets analysis
# gamlssVGD() is a combination of the original TGD()  TGD1() TGD2() functions
# it fits a gamlss model at the training  data (rand==1) and calculates the global
# deviance for the validation set (rand==2)
# TO DO: a) It should also save the residuals for the validation test 
#         Now done the validated residuals are saved as "residVal"
#        b) gamlssCV() needs parallelization
#        c) need a function for comparing CV models    
################################################################################
################################################################################
# functions
################################################################################
#   i) gamlssVGD() for fitting a model and then calculate the deviance for the 
#      extra data
#  ii) VGD() for comparing fitted gamlssVHD models
# iii) getTGD()  after fitting to training data get the global deviance for 
#      the test data
#  iv) TDG() comparing fitted TGD objects
#   v) gamlssCV() for fitting k-folds cross validation
#  vi) CV() comparing fitted CV objects
################################################################################
#rand <- sample(2, 610, replace=T, prob=c(0.6,0.4)
################################################################################
################################################################################
################################################################################
################################################################################
# VALIDATION: this is fitting a gamlss model on a sample from the original data 
#  and calculated the Validated Global Deviance from the new Validated data
gamlssVGD <-function(formula = NULL, 
               sigma.formula = ~1, 
                  nu.formula = ~1, 
                 tau.formula = ~1, 
                        data = NULL, # original data 
                      family = NO,  
                     control = gamlss.control(trace=FALSE),
                        rand = NULL, # 1 for training 2 for validation
                     newdata = NULL, 
                          ...)
 {
################################################################################
# this is to replicate rqres within gamlssVGD enviroment
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
################################################################################
## main function starts here
################################################################################
   sys_call <- sys.call() 
if (is.null(data))   stop("data should be set here")
if (is.null(rand)&&is.null(newdata)) stop("rand or newdata should be set")
if (!is.null(rand))
  {
  if ( any(!rand%in%c(1,2))) stop("rand values should be 1 or 2")
  dataTraining <- subset(data, rand==1)
     dataValid <- subset(data, rand==2) 
}  
       fname <- as.gamlss.family(family)
        dfun <- paste("d", fname$family[[1]],sep="")
        pfun <- paste("p", fname$family[[1]],sep="")
        lpar <- length(fname$parameters)
       dtype <- fname$type
     #  fname$rqres[[1]][["ymin"]]
        if (is.null(formula)) stop("no formula is set in gamlssVGD")        
        if (is.null(data)) stop("the data argument is needed in gamlssVGD")
# FIT MODEL + predict ----------------------------------------------------------
if (!is.null(rand))#  if `rand' is set do this ---------------------------------
 {
          m1 <- gamlss(formula=formula, sigma.formula = sigma.formula, 
                       nu.formula = nu.formula, tau.formula = tau.formula, 
                       data= dataTraining, family=family, control=control, ...) 
     nfitted <- predictAll(m1, newdata=dataValid, data=dataTraining)
dim1newdata  <- dim(dataValid)[1]
 } else
 {                   #if `newdata'  is set do this -----------------------------
          m1 <- gamlss(formula = formula, sigma.formula = sigma.formula, 
                      nu.formula = nu.formula, tau.formula = tau.formula, 
                      data = data, family = family, control = control, ...)
     nfitted <- predictAll(m1, newdata=newdata, data=data)
dim1newdata  <- dim(newdata)[1]
 }
#------------------------------------------------------------------------------
# get  y for new data  if binomial       
if (fname$family[1] %in% .gamlss.bi.list)# if binomial
   {
     if (NCOL(nfitted$y) == 1) 
     {
       y1 <- nfitted$y 
       bd <- nfitted$bd
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
#-------------------------------------------------------------------------------
# jump depending on the number of parameters       
if(lpar==1) 
   {
     if (fname$family[[1]] %in% .gamlss.bi.list)
     {
       devi <-  call(dfun, x=y1, mu =  nfitted$mu, bd=bd, log=TRUE) 
       ures <-  call("rqres", pfun=pfun, type=dtype,
                     ymin=fname$rqres[[1]][["ymin"]], 
                     y=y1, bd=bd, mu= nfitted$mu) 
     } else
     {
       devi <-  call(dfun, x=y1, mu =  nfitted$mu,  log=TRUE) 
       ures <-  call("rqres", pfun=pfun, type=dtype, 
                     ymin=fname$rqres[[1]][["ymin"]], 
                     y=y1, mu= nfitted$mu)
     } 
   }
else if(lpar==2)
   { 
     if (fname$family[[1]] %in% .gamlss.bi.list)
     {
       devi <-  call(dfun, x=y1, mu =  nfitted$mu, sigma= nfitted$sigma,  
                     bd=bd, log=TRUE) 
       ures <-  call("rqres", pfun=pfun, type=dtype, 
                     ymin=fname$rqres[[1]][["ymin"]], 
                     y=y1, mu= nfitted$mu, sigma= nfitted$sigma, bd=bd)
     } else
     {
       devi <-  call(dfun, x=y1, mu =  nfitted$mu, sigma =  nfitted$sigma, log=TRUE) 
       ures <- call("rqres", pfun=pfun,  type=dtype, 
                     ymin=fname$rqres[[1]][["ymin"]], 
                     y=y1, mu= nfitted$mu, sigma= nfitted$sigma )
     } 
   }
else if(lpar==3)
   {
     if (fname$family[[1]] %in% .gamlss.bi.list)
     {
       devi <- call(dfun, x=y1, mu =  nfitted$mu, sigma =  nfitted$sigma, 
                     nu =  nfitted$nu, bd=bd, log=TRUE)
       ures <- call("rqres", pfun=pfun,  type=dtype, 
                    ymin=fname$rqres[[1]][["ymin"]], 
                    y=y1, mu= nfitted$mu, sigma= nfitted$sigma,
                    nu =  nfitted$nu, bd=bd)
     } else
     {
       devi <-  call(dfun, x=y1, mu =  nfitted$mu, sigma =  nfitted$sigma, 
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
       devi <-  call(dfun, x=y1, mu =  nfitted$mu, sigma =  nfitted$sigma, 
                     nu =  nfitted$nu, tau= nfitted$tau, bd=bd, log=TRUE)
       ures <- call("rqres", pfun=pfun,  type=dtype, 
                    ymin=fname$rqres[[1]][["ymin"]], 
                    y=y1, mu= nfitted$mu, sigma= nfitted$sigma,
                    nu =  nfitted$nu, tau = nfitted$tau, bd=bd)
     } else
     {
       devi <-  call(dfun, x=y1, mu =  nfitted$mu, sigma =  nfitted$sigma,
                     nu =  nfitted$nu, tau =  nfitted$tau, log=TRUE)
       ures <- call("rqres", pfun=pfun,  type=dtype, 
                    ymin=fname$rqres[[1]][["ymin"]], 
                    y=y1, mu= nfitted$mu, sigma= nfitted$sigma,
                    nu =  nfitted$nu, tau= nfitted$tau)
      # ures <-  call(pfun, q=y1, mu =  nfitted$mu, sigma =  nfitted$sigma,nu =  nfitted$nu,tau =  nfitted$tau)
     } 
   }
         Vresid <- eval(ures)
       dev.incr <- -2 * eval(devi)
            dev <- sum(dev.incr)
         m1$VGD <- dev
     m1$IncrVGD <- dev.incr
m1$predictError <- dev/dim1newdata
    m1$residVal <- Vresid 
 m1$call$family <- sys_call$family
#   }#----------END of newdata -------------------------------------------------
 class(m1) <- c("gamlssVGD", "gamlss", "gam",    "glm",    "lm"  )
     m1
 } # end of function
################################################################################
################################################################################
################################################################################
################################################################################
VGD <- function(object,...) #UseMethod("AIC")
{
  # local function
is.gamlssVGD <-  function (x)   inherits(x, "gamlssVGD")
  if (length(list(...))) 
  {
      object <- list(object, ...)
    isgamlss <- unlist(lapply(object, is.gamlssVGD))
    if (!any(isgamlss)) stop("some of the objects are not gamlssVGD")
         VGD <- as.numeric(lapply(object, function(x) x$VGD)) 
         val <- cbind(VGD)
        Call <- match.call()
   row.names <- as.character(Call[-1])
       o.val <- order(val)
         val <-  as.data.frame(val[o.val])
       o.r.n <-row.names[o.val]
rownames(val) <- o.r.n 
colnames(val) <- "Pred.GD"
  val
  }
  else 
  { val <- if (is.gamlssVGD(object)) object$VGD 
    else stop(paste("this is not a gamlssVGD object"))
    val 
  }
}
################################################################################
################################################################################
################################################################################
################################################################################   
# this function is modified to accepts gamlss2 objects
# this requires a fitted model and the new data only
getTGD<- function (object,   newdata = NULL, ...) 
{
if (!inherits(object, c("gamlss", "gamlss2"))) stop("not a gamlss object")
if (is.null(newdata)) stop("no newdata is set in VGD")  
if (inherits(object,"gamlss"))
 {
   fname <- as.gamlss.family(object$family[[1]])
    dfun <- paste("d", fname$family[[1]], sep = "")
    pfun <- paste("p", fname$family[[1]],sep="")
    lpar <- length(fname$parameters) 
 nfitted <- predictAll(object, newdata=newdata, ...)
 if (is.null(nfitted$y)) 
   stop("the response variables is missing in the newdata")
 if (fname$family[1] %in% .gamlss.bi.list)# if binomial
   {
     if (NCOL(nfitted$y) == 1) 
     {
       y1 <- nfitted$y 
       bd <- nfitted$bd
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
 if(lpar==1) 
   {
     if (fname$family[[1]] %in% .gamlss.bi.list)
     {
       devi <-  call(dfun, x=y1, mu =  nfitted$mu, bd=bd, log=TRUE) 
       ures <-  call(pfun, q=y1, mu =  nfitted$mu, bd=bd)  
     } else
     {
       devi <-  call(dfun, x=y1, mu =  nfitted$mu,  log=TRUE) 
       ures <-  call(pfun, q=y1, mu =  nfitted$mu)   
     } 
   }
   else if(lpar==2)
   {
     if (fname$family[[1]] %in% .gamlss.bi.list)
     {
       devi <-  call(dfun, x=y1, mu =  nfitted$mu, sigma= nfitted$sigma,  bd=bd, 
                     log=TRUE) 
       ures <-  call(pfun, q=y1, mu =  nfitted$nmu, sigma= nfitted$sigma, bd=bd)  
     } else
     {
       devi <-  call(dfun, x=y1, mu =  nfitted$mu, sigma =  nfitted$sigma, 
                     log=TRUE) 
       ures <-  call(pfun, q=y1, mu =  nfitted$mu, sigma =  nfitted$sigma) 
     } 
   }
   else if(lpar==3)
   {
     if (fname$family[[1]] %in% .gamlss.bi.list)
     {
       devi <-  call(dfun, x=y1, mu =  nfitted$mu, sigma =  nfitted$sigma, nu =  nfitted$nu, bd=bd, log=TRUE)
       ures <-  call(pfun, q=y1, mu =  nfitted$mu, sigma =  nfitted$sigma, nu =  nfitted$nu, bd=bd) 
     } else
     {
       devi <-  call(dfun, x=y1, mu =  nfitted$mu, sigma =  nfitted$sigma, nu =  nfitted$nu, log=TRUE)
       ures <-  call(pfun, q=y1, mu =  nfitted$mu, sigma =  nfitted$sigma, nu =  nfitted$nu)
     } 
   }
   else 
   {
     if (fname$family[[1]] %in% .gamlss.bi.list)
     {
       devi <-  call(dfun, x=y1, mu =  nfitted$mu, sigma =  nfitted$sigma, nu =  nfitted$nu, tau= nfitted$tau, bd=bd, log=TRUE)
       ures <-  call(pfun, q=y1, mu =  nfitted$mu, sigma =  nfitted$sigma, nu =  nfitted$nu, tau= nfitted$tau, bd=bd)
     } else
     {
       devi <-  call(dfun, x=y1, mu =  nfitted$mu, sigma =  nfitted$sigma,nu =  nfitted$nu,tau =  nfitted$tau, log=TRUE)
       ures <-  call(pfun, q=y1, mu =  nfitted$mu, sigma =  nfitted$sigma,nu =  nfitted$nu,tau =  nfitted$tau)
     } 
   }
   Vresid <- qNO(eval(ures))
      dev <- -2*sum(eval(devi))
      out <-list()
      out <- list(TGD=dev,  predictError=dev/dim(newdata)[1], resid=Vresid)   
 } else # if gamlss2
 {
################################################################################   
  get_family <- function (model) # local function
   {
     if  (!missing(model)&&(!inherits(model, c("gamlss", "gamlss2")))) 
       stop("the model should be an gamlss object")  
     if (inherits(model, "gamlss"))
     {
       family <-  if(is.null(model$call$family)) as.gamlss.family(NO) 
                   else as.gamlss.family(model$call$family)
        fname <- model$family[1]  
         type <- family$type
        param <- family$param
        nopar <- family$nopar
         dfun <- paste("d",fname,sep="")
         pfun <- paste("p",fname,sep="")
         qfun <- paste("q",fname,sep="")
        p_d_f <- eval(parse(text=dfun)) 
        c_d_f <- eval(parse(text=pfun))
        q_fun <- eval(parse(text=qfun))
     } else 
     {
       family <- model$family
        fname <- family$family
         type <- family$type
        param <- family$names
        nopar <- length(param)
         dfun <- paste("d",fname,sep="")
         pfun <- paste("p",fname,sep="")
        p_d_f <- eval(family$d) 
        c_d_f <- eval(family$p) 
        q_fun <- eval(family$q) 
     } 
    list(fname=fname, type=type, nopar=nopar, param=param, dfun=dfun, pfun=pfun,
          p_d_f=p_d_f, c_d_f=c_d_f, q_fun=q_fun)  
  }
################################################################################  
  # HERE IS WHAT WE WANT
  # predict in new data and calculate the deviance 
#         param <- as.data.frame(gamlss2:::predict.gamlss2(object, newdata=newdata, type="parameter"))
#          lpar <- length(object$family$names)
#           fam <- get_family(object)
# response.name <- paste(object$call$formula[[2]])
#             y <- newdata[,response.name]
#    logLikIncr <- fam$p_d_f(y, param, log=TRUE)
#           dev <- sum(-2*logLikIncr)
      dev <- -2*logLik(object, newdata=newdata)           
      out <- list(TGD=dev,  predictError=dev/dim(newdata)[1])   
 }   
  class(out) <- "gamlssTGD"
  out
}
################################################################################
################################################################################
################################################################################
################################################################################ 
TGD <- function(object,...) #UseMethod("AIC")
{
# local function
is.TGD <-  function (x)   inherits(x, "gamlssTGD")
if (length(list(...))) 
  {
      object <- list(object, ...)
       isTGD <- unlist(lapply(object, is.TGD))
    if (!any(isTGD)) stop("some of the objects are not TGD")
          TGD <- as.numeric(lapply(object, function(x) x$TGD)) 
          val <- cbind(TGD)
         Call <- match.call()
    row.names <- as.character(Call[-1])
        o.val <- order(val)
         val  <- as.data.frame(val[o.val])
        o.r.n <- row.names[o.val]
rownames(val) <- o.r.n
colnames(val) <- "Pred.GD"
    val
  }
  else 
  {       val <- if (is.TGD(object)) object$TGD 
                  else stop(paste("this is not a TGD object"))
    val 
  }
}
################################################################################
################################################################################
################################################################################
################################################################################
gamlssCV <- function(formula = NULL, 
               sigma.formula = ~1, 
                  nu.formula = ~1, 
                 tau.formula = ~1, 
                        data = NULL, # original data 
                      family = NO,  
                     control = gamlss.control(trace=FALSE),
                      K.fold = 10,
                    set.seed = 123,
                        rand = NULL,
                    parallel = c("no", "multicore", "snow"), 
                       ncpus = 1L, 
                          cl = NULL, ...) 
{
#--------------- PARALLEL-------------------------------------------------------
#----------------SET UP PART----------------------------------------------------
if (missing(parallel)) 
    parallel <- "no"
    parallel <- match.arg(parallel)
     have_mc <- have_snow <- FALSE
if (parallel != "no" && ncpus > 1L) 
{
  if (parallel == "multicore") 
      have_mc <- .Platform$OS.type != "windows"
  else if (parallel == "snow") 
    have_snow <- TRUE
  if (!have_mc && !have_snow) 
        ncpus <- 1L
  loadNamespace("parallel")
}
# -------------- finish parallel------------------------------------------------
#-------------------------------------------------------------------------------
sys_call <- sys.call() 
if (is.null(data))   stop("data should be set here")
       N <- dim(data)[1]
# RNAMES <- rownames(data)
set.seed <- set.seed
    rand <- if (is.null(rand)) sample(K.fold , N, replace=TRUE)
            else rand
if (length(rand)!=N) stop("the length of the rand should be equal to data")
  K.fold <- length(unique(rand))
      CV <- rep(0, K.fold)
 residCV <- rep(0, N)
       i <- sort(unique(rand))
#---------------------------------------
fn <- function(i,...)
 {
    cat("fold ", i, "\n",sep="")
    learn <- gamlssVGD(formula=formula, sigma.formula=sigma.formula, nu.formula=nu.formula, tau.formula=tau.formula,  family=family, control=control, data=data[rand!=i,], newdata=data[rand==i,],...)
    # residCV[rand==i] <<-  learn$residVal
    list(CV=learn$VGD, resid=learn$residVal) 
 }
#pp<-vapply(i, fn, list("gc", "res"))
# ll<-lapply(i, fn)
# ss<-sapply(i, fn)
#========================================
# --------  parallel -----------------------------------------------------------
CV <- if (ncpus > 1L && (have_mc || have_snow)) 
{ 
  if (have_mc) 
  {# sapply(scope, fn)
      parallel::mclapply(i, fn, mc.cores = ncpus)
  }
  else if (have_snow) 
  {
      list(...)
   if (is.null(cl)) 
    {
# make the cluster
# cl <- parallel::makePSOCKcluster(rep("localhost", ncpus))
  cl <- parallel::makeForkCluster(ncpus)
    if (RNGkind()[1L] == "L'Ecuyer-CMRG") 
          parallel::clusterSetRNGStream(cl)
 res <- parallel::parLapply(cl, i, fn)
    parallel::stopCluster(cl)
  res
     }
    else t(parallel::parLapply(cl, i, fn))
    }
} # end parallel ---------------------------------------------------------------
else lapply(i, fn)
  cv <- rep(0, K.fold)
for (i in 1:K.fold)
  {
    residCV[rand==i] <-  CV[[i]]$resid
    cv[i] <- CV[[i]]$CV
  }
 # CV <-  sapply(i, fn) 
  out <- list(CV=sum(cv), allCV=cv, residCV=residCV) 
  class(out) <- "gamlssCV"
  out
}
#-------------------------------------------------------------------------------
CV <- function(object,...) #UseMethod("AIC")
{
  # local function
  is.CV <-  function (x)   inherits(x, "gamlssCV")
  if (length(list(...))) 
  {
    object <- list(object, ...)
    isCV <- unlist(lapply(object, is.CV))
    if (!any(isCV)) stop("some of the objects are not gamlssCV")
    CV <- as.numeric(lapply(object, function(x) x$CV)) 
    val <- cbind(CV)
    Call <- match.call()
    row.names <- as.character(Call[-1])
    o.val <- order(val)
    val  <-  as.data.frame(val[o.val])
    o.r.n <- row.names[o.val]
    rownames(val) <- o.r.n 
    val
  }
  else 
  { val <- if (is.CV(object)) object$CV 
    else stop(paste("this is not a gamlssCV object"))
    val 
  }
}
################################################################################
################################################################################
################################################################################
################################################################################