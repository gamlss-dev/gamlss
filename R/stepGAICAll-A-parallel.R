# trying to make the function parallel
################################################################################
stepGAICAll.A <- function(object, scope = NULL, 
                             sigma.scope = NULL, 
                                nu.scope = NULL, 
                               tau.scope = NULL, 
                                  mu.try = TRUE, 
                               sigma.try = TRUE, 
                                  nu.try = TRUE, 
                                 tau.try = TRUE, 
                               direction = NULL,
                                parallel = c("no", "multicore", "snow"),
                                   ncpus = 1L, 
                                     cl = NULL, 
                          ...)
{
################################################################################
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
if (have_snow)
{
       cl <- parallel::makeForkCluster(ncpus)
  if (RNGkind()[1L] == "L'Ecuyer-CMRG") 
    parallel::clusterSetRNGStream(cl)
  on.exit(parallel::stopCluster(cl))
}         
# -------------- finish parallel------------------------------------------------
################################################################################
################################################################################
#-------------------------------------------------------------------------------
## make sure that the object is visible 
      objectAll <- object
      All.anova <- list()
## if defferent scope is required
    sigma.scope <- if (is.null(sigma.scope)) scope else sigma.scope
       nu.scope <- if (is.null(   nu.scope)) scope else    nu.scope
      tau.scope <- if (is.null(  tau.scope)) scope else   tau.scope  
## get the number of parameters from the distribution
            pdf <- object$family[1]
           npar <- eval(parse(text=pdf))()$nopar
      direction <- if (is.null(direction))
    {
      switch(npar, "both",
             c("forward","both", "backward"),
             c("forward", "forward", "both", "backward", "backward"),
             c("forward", "forward", "forward","both", "backward", "backward", "backward"))
    } else       direction    
if  (npar==1)
{
if (trace)
{
  cat("---------------------------------------------------", "\n") 
}
# get the mu mode
  if ("mu" %in% object$par && mu.try==TRUE)
  {
    current.par <- "mu"
    iferror <- try( assign("objectAll", stepGAIC(objectAll, scope=scope, direction=direction[1], what = "mu",  parallel = parallel,  ncpus = ncpus, cl = cl, ...)), silent = TRUE) 
    if (any(class(iferror)%in%"try-error"))
    { 
      cat("---------------------------------------------------", "\n")
      cat(paste("ERROR: stepGAICAll has failed trying to fit the model for mu forward.", "\n",
                "Maybe the mu model is too complicated for the data.", "\n",
                "The model given is the final before the crush. \n"))
      cat("---------------------------------------------------", "\n")    
      return(objectAll)
    }
    # saving the anova
    All.anova <-list(mu.anova.for=objectAll$anova)
  }
cat("---------------------------------------------------", "\n") 
} 
if  (npar==2)
{
cat("---------------------------------------------------", "\n")
  if ("mu" %in% object$par && mu.try==TRUE)
  {
    current.par <- "mu"
    iferror <- try( assign("objectAll", stepGAIC(objectAll, scope=scope, direction=direction[1], what = "mu",  parallel = parallel,  ncpus = ncpus, cl = cl, ...)), silent = TRUE) 
    if (any(class(iferror)%in%"try-error"))
    { 
      cat("---------------------------------------------------", "\n")
      cat(paste("ERROR: stepGAICAll has failed trying to fit the model for mu forward.", "\n",
                "Maybe the mu model is too complicated for the data.", "\n",
                "The model given is the final before the crush. \n"))
      cat("---------------------------------------------------", "\n")    
      return(objectAll)
    }
# saving the anova
    All.anova <-list(mu.anova.for=objectAll$anova)
  }
  if ("sigma" %in% object$par && sigma.try==TRUE)
  {
    cat("---------------------------------------------------", "\n")
    current.par <- "sigma"
    iferror <- try( assign("objectAll", stepGAIC(objectAll, scope=sigma.scope, direction=direction[2], what = "sigma", parallel = parallel,  ncpus = ncpus, cl = cl, ...))  , silent = TRUE) 
    if (any(class(iferror)%in%"try-error"))
    {
      cat("---------------------------------------------------", "\n") 
      cat(paste("ERROR: stepGAICAll has failed trying to fit the model for sigma forward.", "\n",
                "Maybe the sigma model is too complicated for the data.", "\n",
                "The model given is the final before the crush. \n"))    
      cat("---------------------------------------------------", "\n")
      return(objectAll) 
    } 
    All.anova$sigma.anova <- objectAll$anova      
  }
  if ("mu" %in% object$par && mu.try==TRUE && current.par!="mu")
  {
    cat("---------------------------------------------------", "\n")  
    current.par <- "mu"
    iferror <- try( assign("objectAll", stepGAIC(objectAll, scope=scope, direction=direction[3],  what = "mu", parallel = parallel,  ncpus = ncpus, cl = cl, ...))  , silent = TRUE) 
    if (any(class(iferror)%in%"try-error"))
    { 
      cat("---------------------------------------------------", "\n") 
      cat(paste("ERROR: stepGAICAll has failed trying to fit the model for mu backward.", "\n",
                "Maybe the mu model is too complicated for the data,", "\n",
                "The model given is the final before the crush. \n"))
      cat("---------------------------------------------------", "\n")    
      return(objectAll) 
    } 
    All.anova$mu.anova.back<-objectAll$anova        
}
  cat("---------------------------------------------------", "\n")    
  objectAll$anovaAll<-All.anova    
} 
if (npar==3)
{
  cat("---------------------------------------------------", "\n")
  if ("mu" %in% object$par && mu.try==TRUE)
  {
    current.par <- "mu"
    iferror <- try( assign("objectAll", stepGAIC(objectAll, scope=scope, direction=direction[1], what = "mu",  parallel = parallel,  ncpus = ncpus, cl = cl, ...)), silent = TRUE) 
    if (any(class(iferror)%in%"try-error"))
    { 
      cat("---------------------------------------------------", "\n")
      cat(paste("ERROR: stepGAICAll has failed trying to fit the model for mu forward.", "\n",
                "Maybe the mu model is too complicated for the data.", "\n",
                "The model given is the final before the crush. \n"))
      cat("---------------------------------------------------", "\n")    
      return(objectAll)
    }
    # saving the anova
    All.anova <-list(mu.anova.for=objectAll$anova)
  }
  if ("sigma" %in% object$par && sigma.try==TRUE)
  {
    cat("---------------------------------------------------", "\n")
    current.par <- "sigma"
    iferror <- try( assign("objectAll", stepGAIC(objectAll, scope=sigma.scope, direction=direction[2], what = "sigma", parallel = parallel,  ncpus = ncpus, cl = cl, ...))  , silent = TRUE) 
    if (any(class(iferror)%in%"try-error"))
    {
      cat("---------------------------------------------------", "\n") 
      cat(paste("ERROR: stepGAICAll has failed trying to fit the model for sigma forward.", "\n",
                "Maybe the sigma model is too complicated for the data.", "\n",
                "The model given is the final before the crush. \n"))    
      cat("---------------------------------------------------", "\n")
      return(objectAll) 
    } 
    All.anova$sigma.anova.for <- objectAll$anova      
  }
  if (  "nu" %in% object$par && nu.try==TRUE)
  {
    cat("---------------------------------------------------", "\n") 
    current.par <- "nu"
    iferror <- try( assign("objectAll", obj <- stepGAIC(objectAll, scope=nu.scope,direction=direction[3], what = "nu",parallel = parallel,  ncpus = ncpus, cl = cl, ...)), silent = TRUE)           
    if (any(class(iferror)%in%"try-error"))
    { 
      cat("---------------------------------------------------", "\n")
      cat(paste("ERROR: stepGAICAll has failed trying to fit the model for nu forward.", "\n",
                "Maybe the nu model is too complicated for the data.", "\n",
                "The model given is the final before the crush. \n"))    
      cat("---------------------------------------------------", "\n")
      return(objectAll) 
    }  
    All.anova$nu.anova.for <- objectAll$anova      
  } 
}  
if (npar==4)
{
cat("---------------------------------------------------", "\n")
  if ("mu" %in% object$par && mu.try==TRUE)
  {
    current.par <- "mu"
    iferror <- try( assign("objectAll", stepGAIC(objectAll, scope=scope, direction=direction[1], what = "mu",  parallel = parallel,  ncpus = ncpus, cl = cl, ...)), silent = TRUE) 
    if (any(class(iferror)%in%"try-error"))
    { 
      cat("---------------------------------------------------", "\n")
      cat(paste("ERROR: stepGAICAll has failed trying to fit the model for mu forward.", "\n",
                "Maybe the mu model is too complicated for the data.", "\n",
                "The model given is the final before the crush. \n"))
      cat("---------------------------------------------------", "\n")    
      return(objectAll)
    }
    # saving the anova
    All.anova <-list(mu.anova.for=objectAll$anova)
  }
#eval(expression(objectAll$PPP<-objectAll$anova), envir=.GlobalEnv)
# get the sigma model
#-------------------- 
  if ("sigma" %in% object$par && sigma.try==TRUE)
  {
    cat("---------------------------------------------------", "\n")
    current.par <- "sigma"
    iferror <- try( assign("objectAll", stepGAIC(objectAll, scope=sigma.scope, direction=direction[2], what = "sigma", parallel = parallel,  ncpus = ncpus, cl = cl, ...))  , silent = TRUE) 
    if (any(class(iferror)%in%"try-error"))
    {
      cat("---------------------------------------------------", "\n") 
      cat(paste("ERROR: stepGAICAll has failed trying to fit the model for sigma forward.", "\n",
                "Maybe the sigma model is too complicated for the data.", "\n",
                "The model given is the final before the crush. \n"))    
      cat("---------------------------------------------------", "\n")
      return(objectAll) 
    } 
    All.anova$sigma.anova.for <- objectAll$anova      
  }
# get the nu model
  #------------------ 
  if (  "nu" %in% object$par && nu.try==TRUE)
  {
    cat("---------------------------------------------------", "\n") 
    current.par <- "nu"
    iferror <- try( assign("objectAll", obj <- stepGAIC(objectAll, scope=nu.scope,direction=direction[3], what = "nu",parallel = parallel,  ncpus = ncpus, cl = cl, ...)), silent = TRUE)           
    if (any(class(iferror)%in%"try-error"))
    { 
      cat("---------------------------------------------------", "\n")
      cat(paste("ERROR: stepGAICAll has failed trying to fit the model for nu forward.", "\n",
                "Maybe the nu model is too complicated for the data.", "\n",
                "The model given is the final before the crush. \n"))    
      cat("---------------------------------------------------", "\n")
      return(objectAll) 
    }  
    All.anova$nu.anova.for <- objectAll$anova      
  }
  # get the tau model
  #------------------           
  if ( "tau" %in% object$par && tau.try==TRUE)
  {
    cat("---------------------------------------------------", "\n") 
    current.par <- "tau"
    iferror <- try( assign("objectAll", stepGAIC(objectAll, scope=tau.scope,direction=direction[4], what = "tau", parallel = parallel,  ncpus = ncpus, cl = cl, ...)), silent = TRUE) 
    if (any(class(iferror)%in%"try-error"))
    { 
      cat("---------------------------------------------------", "\n")
      cat(paste("ERROR: stepGAICAll has failed trying to fit the model for tau forward.", "\n",
                "Maybe the tau model is too complicated for the data.", "\n",
                "The model given is the final before the crush. \n"))    
      cat("---------------------------------------------------", "\n")
      return(objectAll) 
    }
    All.anova$tau.anova.for<-objectAll$anova       
  }   
  # get the nu model
  #------------------ 
  if (  "nu" %in% object$par && nu.try==TRUE && current.par!="nu")
  { 
    cat("---------------------------------------------------", "\n") 
    current.par <- "nu"
    iferror <- try(assign("objectAll", stepGAIC(objectAll, scope=nu.scope, direction=direction[5], what = "nu", parallel = parallel,  ncpus = ncpus, cl = cl, ...))  , silent = TRUE) 
    if (any(class(iferror)%in%"try-error"))
    { 
      cat("---------------------------------------------------", "\n")
      cat(paste("ERROR: stepGAICAll has failed trying to fit the model for nu backward.", "\n",
                "Maybe the nu model is too complicated for the data.", "\n",
                "The model given is the final before the crush. \n"))    
      cat("---------------------------------------------------", "\n")
      return(objectAll) 
    }  
    All.anova$nu.anova.back<-objectAll$anova
  }
  # get the sigma model
  #------------------ 
  if ("sigma" %in% object$par && sigma.try==TRUE && current.par!="sigma")
  {
    cat("---------------------------------------------------", "\n")
    current.par <- "sigma"
    iferror <- try( assign("objectAll", stepGAIC(objectAll, scope=sigma.scope, direction=direction[6], what = "sigma", parallel = parallel,  ncpus = ncpus, cl = cl, ...)) , silent = TRUE) 
    if (any(class(iferror)%in%"try-error"))
    { cat("---------------------------------------------------", "\n")
      cat(paste("ERROR: stepGAICAll has failed trying to fit the model for sigma backward.", "\n",
                "Maybe the sigma model is too complicated for the data.", "\n",
                "The model given is the final before the crush. \n"))    
      cat("---------------------------------------------------", "\n")
      return(objectAll) 
    }  
    All.anova$sigma.anova.back <- objectAll$anova
  }
  # get the final mu model
  #------------------ 
  if ("mu" %in% object$par && mu.try==TRUE && current.par!="mu")
  {
    cat("---------------------------------------------------", "\n")  
    current.par <- "mu"
    iferror <- try( assign("objectAll", stepGAIC(objectAll, scope=scope, direction=direction[7],  what = "mu", parallel = parallel,  ncpus = ncpus, cl = cl, ...))  , silent = TRUE) 
    if (any(class(iferror)%in%"try-error"))
    { 
      cat("---------------------------------------------------", "\n") 
      cat(paste("ERROR: stepGAICAll has failed trying to fit the model for mu backward.", "\n",
                "Maybe the mu model is too complicated for the data,", "\n",
                "The model given is the final before the crush. \n"))
      cat("---------------------------------------------------", "\n")    
      return(objectAll) 
    } 
    All.anova$mu.anova.back<-objectAll$anova        
  }
  cat("---------------------------------------------------", "\n")    
  objectAll$anovaAll<-All.anova  
}

objectAll
}
