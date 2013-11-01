stepGAICAll.A <- function(object, scope=NULL, sigma.scope=NULL, nu.scope=NULL, tau.scope=NULL,
                                  mu.try=TRUE, sigma.try=TRUE, nu.try=TRUE, tau.try=TRUE, ...)
{
    ## make sure that the object is visible 
     objectAll<- object
    # on.exit(rm(objectAll, envir=.GlobalEnv)) #' delete on exit
#    env <- attach(NULL, name="Object_All_Env")
#    assign("objectAll", object, envir=env)
#    on.exit(detach(Object_All_Env))
    All.anova <-list()
    ## if defferent scope is required
    sigma.scope <- if (is.null(sigma.scope)) scope else sigma.scope
       nu.scope <- if (is.null(   nu.scope)) scope else    nu.scope
      tau.scope <- if (is.null(  tau.scope)) scope else   tau.scope  
    # get the mu model
    #-----------------
    cat("---------------------------------------------------", "\n")
     if ("mu" %in% object$par && mu.try==TRUE)
         {
     current.par <- "mu"
         iferror <- try( assign("objectAll", stepGAIC.VR(objectAll, scope=scope, direction="forward", what = "mu", ...)), silent = TRUE) 
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
        iferror <- try( assign("objectAll", stepGAIC.VR(objectAll, scope=sigma.scope, direction="forward", 
                       what = "sigma", ...))  , silent = TRUE) 
             if (any(class(iferror)%in%"try-error"))
                 {
                  cat("---------------------------------------------------", "\n") 
                  cat(paste("ERROR: stepGAICAll has failed trying to fit the model for sigma forward.", "\n",
                     "Maybe the sigma model is too complicated for the data.", "\n",
                     "The model given is the final before the crush. \n"))    
                  cat("---------------------------------------------------", "\n")
                   return(objectAll) 
                 } 
            All.anova$sigma.anova.for<-objectAll$anova      
        }
     # get the nu model
     #------------------ 
     if (  "nu" %in% object$par && nu.try==TRUE)
        {
         cat("---------------------------------------------------", "\n") 
          current.par <- "nu"
        iferror <- try( assign("objectAll", obj <- stepGAIC.VR(objectAll, scope=nu.scope, direction="forward", 
                         what = "nu", ...)), silent = TRUE)           
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
        iferror <- try( assign("objectAll", stepGAIC.VR(objectAll, scope=tau.scope, direction="forward", 
                       what = "tau", ...)), silent = TRUE) 
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
        iferror <- try(assign("objectAll", stepGAIC.VR(objectAll, scope=nu.scope, direction="backward", 
                        what = "nu", ...))  , silent = TRUE) 
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
        iferror <- try( assign("objectAll", stepGAIC.VR(objectAll, scope=sigma.scope, direction="backward", 
                        what = "sigma", ...)) , silent = TRUE) 
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
      iferror <- try( assign("objectAll", stepGAIC.VR(objectAll, scope=scope, direction="backward", 
                     what = "mu", ...))  , silent = TRUE) 
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
objectAll
}
