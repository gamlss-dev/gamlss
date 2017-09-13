##-----------------------------------------------------------------------------
#------------------------------------------------------------------------------- 
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# this the predictAll() new function 
# allows the user to get all the parameters using the  predict.gamlss().
# creates a list containing y if exist in the newdata  and predicted values for 
#  mu sigma nu and tau
# # note that default type is "response"
# an option of  new version of predictAll() is  based on weighrs
# The idea is that the fitted model is refitted with zero weights for the 
# predictive valus and that allows to have also stadard errors for prediction.
# refit the whole model could be slow
# the output can be a list or a matrix with mu, sigma, nu and tau
# and se's
# The use.weights=T  function combines the data and refit the model with weights
# to get the predictions 
# ##------------------------------------------------------------------------------
# ## to do
# ##   i) type should  c("link", "response") is working for use.weights=T
#            but I am not sure about term
# ##  ii) if explanatory variables are not in the list 
# ##        we should use the means or for factors the most common level (not implemented)
# ## iii) we should allow lists and data.frame (rather than only data.frame) 
# ##       Probably not
# ##  iv) chech for exceptions for binary and survival data 
# ##  iiv) need a option of how to get the y-variable prediction it is crusial 
# ##       when estimating the deviance in the fit OK see y. 
#----------------------------------------------------------------------------------------
predictAll <-function(object, 
                    newdata = NULL, 
                       type = c("response", "link", "terms"),# note that default is "response" 
                      terms = NULL,   
                     se.fit = FALSE, 
                use.weights = FALSE, 
                       data = NULL, 
                    y.value = "median", # a function or numeric values
                     set.to = .Machine$double.xmin,# .Machine$double.eps^.8,
                     output = c("list", "matrix"), # how to save it
                      ...)       
{
##------------------------------------------------------------------------------
##------------------------------------------------------------------------------
##-------- concat starts here---------------------------------------------------
  concat <- function(..., names=NULL) 
  { 
    tmp <- list(...) 
    if(is.null(names)) names <- names(tmp) 
    if(is.null(names)) names <- sapply( as.list(match.call()), deparse)[-1] 
    if( any( 
      sapply(tmp, is.matrix) 
      | 
      sapply(tmp, is.data.frame) ) ) 
    { 
      len <- sapply(tmp, function(x) c(dim(x),1)[1] ) 
      len[is.null(len)] <- 1 
      data <- rbind( ... ) 
    } 
    else 
    { 
      len <- sapply(tmp,length) 
      data <- unlist(tmp)     
    } 
    namelist <- factor(rep(names, len), levels=names)          
    return( data.frame( data, source=namelist) ) 
  } 
##----------concat finish here-------------------------------------------------
##-----------------------------------------------------------------------------
## main function starts here
##-----------------------------------------------------------------------------  
## structure 
## if no new data get the fitted values for all 
##  else if newdata and use.weights==FALSE -> predict individually
##       if newdat  and use.weights == TRUE  -> refit the model     
   type <- match.arg(type)
 output <- match.arg(output)
   Call <- object$call # the call
      N <- object$N    # the number of observations
   npar <- length(object$parameters)
##       geting the  response as character 
 ResCha <- as.character(object$mu.formula[2][[1]])
##
 DatA <- if (is.null(data))
 {        ## if it is not provided then get it from the original call
   if (!is.null(Call$data)) eval(Call$data) 
   else stop("define the original data using the option data") 
 }
 else data # if it provide get it 
## -----------------------------------------------------------------------------
if (is.null(newdata)) # if no new data then give all the fitted from the old predict
{
  if (output=="list")
   {
    out <- list(y=object$y)
    if ("mu" %in% object$par)  
      out$mu <- lpred(object, what = "mu", type = type, terms = terms, se.fit = se.fit, ... )
    if ("sigma" %in% object$par)  
      out$sigma <- lpred(object, what = "sigma", type = type, terms = terms, se.fit = se.fit, ... )
    if (  "nu" %in% object$par)     
      out$nu <- lpred(object, what = "nu", type = type, terms = terms, se.fit = se.fit, ... )
    if ( "tau" %in% object$par)    
      out$tau <- lpred(object, what = "tau", type = type, terms = terms, se.fit = se.fit, ... )
    attr(out, "family") <- object$family
    return(out)
   }
  if (output=="matrix")
   {
            prematrix <- matrix(0, ncol=length(object$parameters)+1, nrow=N)
  colnames(prematrix) <- c("y", object$parameters)
        prematrix[,1] <- object$y
     if (se.fit)
      { 
                 predse <- matrix(0, ncol=length(object$parameters), nrow=N) 
       colnames(predse) <- object$parameters
       for (i in 1:length(object$parameters))
        { 
               PVaSE <-  lpred(object, parameter=object$parameter[i], type=type, se.fit=TRUE)
     prematrix[,i+1] <-  PVaSE$fit  
          predse[,i] <- PVaSE$se.fit  
        }
       colnames(predse) <- paste(colnames(predse), ".se", sep="")
       matout <-  cbind(prematrix, predse) 
       attr(matout, "family") <- object$family
       return(matout)
     } else
     {  for (i in object$parameters)
       { prematrix[,i] <- lpred(object, parameter=i, type=type)  }
       attr(prematrix, "family") <- object$family
       return(prematrix)  
     }    
  }   
}# end  if (is.null(newdata)) 
  else #         --------   if new data -------- 
  {
    pN <- dim(newdata)[1]
if (!(inherits(newdata, "data.frame")))
      stop("newdata must be a data frame ") # or a frame mumber
    
if ((use.weights==FALSE)&&(se.fit==FALSE))#  if use.weights is FALSE and se.fit=FALSE 
{                                 # use predict()
  
  if (output=="list")
  {
    out <- list()
    if ("mu" %in% object$par) #
      out$mu <- predict(object,newdata=newdata, what = "mu", type = type, terms = terms, se.fit = se.fit, data=DatA )
    if ("sigma" %in% object$par)  
      out$sigma <- predict(object, newdata=newdata, data=DatA,   what = "sigma", type = type, terms = terms, se.fit = se.fit)
    if ("nu" %in% object$par)  
      out$nu <- predict(object, newdata=newdata, data=DatA,  what = "nu", type = type, terms = terms, se.fit = se.fit )
    if ("tau" %in% object$par)  
      out$tau <- predict(object, newdata=newdata, data=DatA , what = "tau", type = type, terms = terms, se.fit = se.fit)
    if (as.character(object$mu.formula[[2]])%in%names(newdata)) 
    out$y <-  newdata[,as.character(object$mu.formula[[2]])]
    attr(out, "family") <- object$family
    return(out)    
  }
   if (output=="matrix")
   {
     ifY <- as.character(object$mu.formula[[2]])%in%names(newdata)
     if (ifY)
     {
          prematrix <- matrix(0, ncol=length(object$parameters)+1, nrow=pN)
colnames(prematrix) <- c("y", object$parameters)
     prematrix[,1]  <-  newdata[,as.character(object$mu.formula[[2]])]
     }
     else 
     {
          prematrix <- matrix(0, ncol=length(object$parameters), nrow=pN)
colnames(prematrix) <- object$parameters
     }
     for (i in object$parameters)
     { prematrix[,i] <-  predict(object,newdata=newdata, what = i, type = type, terms = terms, se.fit = se.fit, data=DatA )}
     attr(prematrix, "family") <- object$family 
  return(prematrix) 
   }
}
if ((use.weights==TRUE)||(se.fit==TRUE))
   {
     # DatA <- if (is.null(data))
     # {        ## if it is not provided then get it from the original call
     #   if (!is.null(Call$data)) eval(Call$data) 
     #   else stop("define the original data using the option data") 
     # }
     # else data # if it provide get it 
## here we should check whether newdata has ResCha
       ifY <- ResCha%in%names(newdata)     
##       keep only the same variables 
##       this assumes that all the relevant variables will be in newdata
##       what happens if not?     
##      The question is if newdata does not contain all variables what we do 
     DatA <- DatA[match(names(newdata),names(DatA), nomatch = 0)]# mikis 
## merge the two data together
     DatA <- concat(DatA,newdata)
## get the old weights
     oldweights <-object$weights
## get the old response     
     oldresponse <- object$y
## create the new weights    
     pN <- dim(newdata)[1]
     DatA$Weights <- c( rep(1, N)*oldweights, rep(set.to, pN))
## now ifY is true then yvar should take values from the newdata
  yval <-  if (ifY) newdata[,ResCha] 
          else{ # otherwise use y.value
            if (!is.numeric(y.value))  rep(do.call(y.value, list(x=oldresponse)), pN)
            else rep(y.value, pN)  
          } 
##   cat("y", yval)
     DatA[[ResCha ]] <- c(oldresponse, yval)
## fits the model
     newobj <-  update(object, data=DatA, weights=DatA$Weights, trace=FALSE)
     # needs some checking here   
     DiffAIC  <- AIC(object, newobj, k=0)
     if (abs(diff(DiffAIC[,2]))>1)   warning("the global deviances differ more that 1 ")
     if (abs(diff(DiffAIC[,1]))>.5)  warning("the degrees of freedom differ more that 0.5 ")
     cat("dev:",deviance(object),deviance(newobj), deviance(object)-deviance(newobj),"\n")
     cat("df:",object$df.fit,newobj$df.fit, object$df.fit-newobj$df.fit,"\n")
## preper the output
     prematrix <- matrix(0, ncol=length(object$parameters)+1, nrow=pN)
     colnames(prematrix) <- c("y",  object$parameters)
     prematrix[, "y"] <- yval
     if ( se.fit)
     { 
       predse <- matrix(0, ncol=length(object$parameters), nrow=pN) 
       colnames(predse) <- object$parameters
       for (i in object$parameters)
       { 
         PVaSE <-  lpred(newobj, parameter=i, type=type, se.fit=TRUE)
         prematrix[,i] <-  PVaSE$fit[(N+1):(N+pN)]  
         predse[,i] <- PVaSE$se.fit[(N+1):(N+pN)]  
       }
       colnames(predse) <- paste(colnames(predse), ".se", sep="")
       listout <-  if (output=="matrix") cbind(prematrix, predse) else as.list(data.frame(cbind(prematrix, predse)))
       return(listout)
     } else
     { 
       for (i in object$parameters)
       { prematrix[,i] <- lpred(newobj, parameter=i, type=type)[(N+1):(N+pN)]  }
       listout <- if (output=="matrix") prematrix else as.list(data.frame(prematrix))  
       #  lapply(seq_len(ncol(prematrix)), function(i) prematrix[,i]) 
       return(listout)  
     }    
   }
  }  
}
#------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# ## the new version of predict based on weighrs
# ## The idea is that the fitted model is refitted with zero weights for the 
# ## predictive valus that allows to have also stadard errors for prediction
# ## No paameters values are needed here since we refit the whole model
# ## the output should be a matrix with mu, sigma, nu and tau
# ## or a list with two matrices fit and se's
# ## the function combines the data
# ## refit the model with weights
# ## and get the predictions 
# ##------------------------------------------------------------------------------
# ## to do
# ##    i)  Newdata should be compulsory OK
# ##   ii) type should  c("link", "response") but I am not sure about term
# ##  iii) if explanatory variables are not in the list 
# ##        we should use the means or for factors the most common level 
# ##    v) we should allow lists and data.frame (rather than only data.frame) 
# ##       Probably not
# ##   vi) Exceptions for binary and survival data
# ##  vii) need a option of how to get the y-variable prediction it is crusial 
# ##       when estimating the deviance in the fit 
# ## 
# #-------------------------------------------------------------------------------
# #-------------------------------------------------------------------------------
# predictMe <- function(object, 
#                    newdata = NULL, # data.frame
#                       type = c("response", "link", "terms"), # terms not working 
#                      terms = NULL, 
#                     se.fit = FALSE, 
#                       data = NULL, # you need an option how
#                    y.value = "median", # a function or numeric
#                      set.to = .Machine$double.xmin,# .Machine$double.eps^.8,
#                 if.matrix = FALSE,
#                       ...)
# {
# ## this little function put data frames together 
# ##  originated from an the R-help reply by B. Ripley
# ##------------------------------------------------------------------------------
# ##------------------------------------------------------------------------------
# ##-------- concat starts here---------------------------------------------------
# concat <- function(..., names=NULL) 
#   { 
#     tmp <- list(...) 
#     if(is.null(names)) names <- names(tmp) 
#     if(is.null(names)) names <- sapply( as.list(match.call()), deparse)[-1] 
#     if( any( 
#             sapply(tmp, is.matrix) 
#             | 
#             sapply(tmp, is.data.frame) ) ) 
#       { 
#         len <- sapply(tmp, function(x) c(dim(x),1)[1] ) 
#         len[is.null(len)] <- 1 
#         data <- rbind( ... ) 
#       } 
#     else 
#       { 
#         len <- sapply(tmp,length) 
#         data <- unlist(tmp)     
#       } 
#     namelist <- factor(rep(names, len), levels=names)          
#     return( data.frame( data, source=namelist) ) 
#   } 
# ##----------concat finish here--------------------------------------------------
# ##-----------------------------------------------------------------------------
# ##-----------------------------------------------------------------------------
# ##   main function starts here
# ##------------------------------------------------------------------------------
# ## If no new data stop
# if (is.null(newdata))  # 
#     {
#    stop("newdata is needed for this function")
#     }
# # only data frame 
#   if (!(inherits(newdata, "data.frame")))
#     stop("newdata must be a data frame ") # or a frame mumber
# ## getting  type   
#        type <- match.arg(type)
# ## get the original call 
#        Call <- object$call
# ## we need both the old and the new data sets
# ## the argument data can be provided by predict
# data <- if (is.null(data))
# {        ## if it is not provided then get it from the original call
#             if (!is.null(Call$data)) eval(Call$data) 
#             else stop("define the original data using the option data") 
#            }
#         else data # if it provide get it 
# ## keep only the same variables 
# ## this assumes that all the relevant variables will be in newdata
# ## what happens if not?
#       ResCha <- as.character(object$mu.formula[2][[1]])
# # The question is if newdata does not contain all variables what we do 
#              data <- data[match(names(newdata),names(data))]  
# ## merge the two data together
#             data <- concat(data,newdata)
# ## get the old weights
#      oldweights <-object$weights
#     oldresponse <- object$y
# ## create the new weights    
#               N <- object$N    
#              pN <- dim(newdata)[1]
#    data$Weights <- c( rep(1, N)*oldweights, rep(set.to, pN))
#    yval <-  if (!is.numeric(y.value))  rep(do.call(y.value, list(x=oldresponse)), pN)
#             else rep(y.value, pN)  
# ##   cat("y", yval)
# data[[ResCha ]] <- c(oldresponse, yval)
# # fits the model
#          newobj <-  update(object, data=data, weights=Weights, trace=FALSE)
#  # needs some checking here   
#        DiffAIC  <- AIC(object, newobj, k=0)
#    if (abs(diff(DiffAIC[,2]))>1)   warning("the global deviances differ more that 1 ")
#    if (abs(diff(DiffAIC[,1]))>.5)  warning("the degrees of freedom differ more that 0.5 ")
#        cat("dev:",deviance(object),deviance(newobj), deviance(object)-deviance(newobj),"\n")
#        cat("df:",object$df.fit,newobj$df.fit, object$df.fit-newobj$df.fit,"\n")
# # preper the output     
#                  prematrix <- matrix(0, ncol=length(object$parameters), nrow=pN)
#        colnames(prematrix) <- object$parameters
# if ( se.fit)
# { 
#            predse <- matrix(0, ncol=length(object$parameters), nrow=pN) 
#  colnames(predse) <- object$parameters
#   for (i in object$parameters)
#   { 
#             PVaSE <-  lpred(newobj, parameter=i, type=type, se=TRUE)
#     prematrix[,i] <-  PVaSE$fit[(N+1):(N+pN)]  
#        predse[,i] <- PVaSE$se.fit[(N+1):(N+pN)]  
#   }
#         
#  colnames(predse) <- paste(colnames(predse), ".se", sep="")
#           listout <-  if (if.matrix) cbind(prematrix, predse) else as.list(data.frame(cbind(prematrix, predse)))
#           listout
# } else
# { 
#   for (i in object$parameters)
#   { prematrix[,i] <- lpred(newobj, parameter=i, type=type)[(N+1):(N+pN)]  }
#   listout <- if (if.matrix) prematrix else as.list(data.frame(prematrix))  
#   #  lapply(seq_len(ncol(prematrix)), function(i) prematrix[,i]) 
#   returm(listout)  
# }    
#  }
# #-------------------------------------------------------------------------------    
# #-------------------------------------------------------------------------------  
predictAll2 <-function(object,
                    newdata = NULL,
                       type = c("response", "link", "terms"),# note that default is "response"
                      terms = NULL,
                     se.fit = FALSE,
                             ...)
 {
  type <- match.arg(type)
## if no new data then give all the fitted from the old
 if (is.null(newdata))  #
    {
    out <- list(y=object$y)
      if ("mu" %in% object$par)
         out$mu <- lpred(object, what = "mu", type = type, terms = terms, se.fit = se.fit, ... )
     if ("sigma" %in% object$par)
      out$sigma <- lpred(object, what = "sigma", type = type, terms = terms, se.fit = se.fit, ... )
     if (  "nu" %in% object$par)
         out$nu <- lpred(object, what = "nu", type = type, terms = terms, se.fit = se.fit, ... )
     if ( "tau" %in% object$par)
        out$tau <- lpred(object, what = "tau", type = type, terms = terms, se.fit = se.fit, ... )
        attr(out, "family") <- object$family
    return(out)
    }
  else
    {
      out <- list()
      if ("mu" %in% object$par) #
         out$mu <- predict(object,newdata=newdata, what = "mu", type = type, terms = terms, se.fit = se.fit, ... )
     if ("sigma" %in% object$par)
      out$sigma <- predict(object, newdata=newdata, what = "sigma", type = type, terms = terms, se.fit = se.fit, ... )
     if ("nu" %in% object$par)
         out$nu <- predict(object, newdata=newdata, what = "nu", type = type, terms = terms, se.fit = se.fit, ... )
     if ("tau" %in% object$par)
        out$tau <- predict(object, newdata=newdata, what = "tau", type = type, terms = terms, se.fit = se.fit, ... )
     if (as.character(object$mu.formula[[2]])%in%names(newdata))
          out$y <-  newdata[,as.character(object$mu.formula[[2]])]
     attr(out, "family") <- object$family
     #out<- list(out,  family=object$family, parameters=object$parameters,  call=object$call,
     #         weights=object$weights, G.deviance=object$G.deviance, N=object$N, type=object$type,
     #         #residuals=object ,
     #         noObs=object$noObs,df.fit=object$df.fit, df.residual=object$df.residuals)
     #      class(out) <- c("gamlssPredict", "gamlss")
       return(out)
     }
 }
#----------------------------------------------------------------------------------------

