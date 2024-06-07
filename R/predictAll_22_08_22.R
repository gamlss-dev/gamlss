###############################################################################
###############################################################################
###############################################################################
###############################################################################

# last change 7-6-24
# Tim Cole suggestion line 158-159
# this the predictAll() function 
# allows the user to get all the parameters using the  predict.gamlss().
# creates a list containing y if exist in the newdata and predicted values for 
# mu sigma nu and tau
# Note that default type is "response" not "link'
# an option of new version of predictAll() is based on weights
# The idea is that the fitted model is refitted with zero weights for the 
# predictive valus and that allows to have also standard errors for prediction.
# refit the whole model could be slow
# the output can be a list or a matrix with mu, sigma, nu and tau
# and se's
# The use.weights=T  function combines the data and refit the model with weights
# to get the predictions 
# the binomial error is NOT checked in the new revised method 
# The revised method allows Surv() response variables 
###-------------------------------------------------------------------------
# ## to do
# ##   i) type should  c("link", "response") is working for use.weights=T
#            but I am not sure about term
# ##  ii) if explanatory variables are not in the list 
# ##        we should use the means or for factors the most common level as defaults
#           (not implemented)
# ## iii) we should allow lists and data.frame (rather than only data.frame) 
# ##       Probably not
# ##  iv) check for exceptions for binary and survival data
#          (This important 27-02-21) 
# ##  iiv) need a option of how to get the y-variable prediction it is crucial 
# ##       when estimating the deviance in the fit OK see y. 
#----------------------------------------------------------------------------------------
predictAll <-function(object, 
                    newdata = NULL, 
                       type = c("response", "link", "terms"),
                                      # note that default is "response" 
                      terms = NULL,   
                     se.fit = FALSE, 
                use.weights = FALSE, 
                       data = NULL, 
                    y.value = "median", # a function or numeric values
                     set.to = .Machine$double.xmin,# .Machine$double.eps^.8,
                     output = c("list","data.frame", "matrix"), # how to save it
                      ...)       
{
################################################################################
################################################################################
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
###############################################################################
###############################################################################
###############################################################################
## main function starts here
###############################################################################  
## structure there are three cases
## i) if no new data (with or without se.fit) get the fitted values for each param 
## ii) a) newdata without se.fit get the fitted values for each param 
##     b) newdata and  se.fit (use.weights == TRUE to  refit the model)
##  
############################################################################### 
gamlss.bi.list <- .gamlss.bi.list          
   type <- match.arg(type) # 
 output <- match.arg(output)
   Call <- object$call # the call
      N <- object$N    # the number of observations
  fname <- object$family[1]        
   npar <- length(object$parameters)
## geting the  response as character 
 ResCha <- if(is.null(names(attr(object$mu.formula,"dataClass")[1]))) as.character(object$mu.formula[2][[1]])
       else names(attr(object$mu.formula,"dataClass")[1])
 if (any(ResCha=="Surv")){
   ResCha <- ResCha[-which(ResCha=="Surv")]
 }
    lRes <- length(ResCha)
## data
    DatA <- if (is.null(data))
 {   ## if it is not provided then get it from the original call
   if (!is.null(Call$data)) eval(Call$data) 
   else stop("define the original data using the option data") 
 }
 else data # if it provide get it 
############################################################################
## if no new data then give all the fitted from the old predict
## and ger also se.fit if true
############################################################################ 
if (is.null(newdata)) # OLD DATA (CASE (i))
{ # OLD DATA STARS
  out <- list(y=DatA[,ResCha])
    if ("mu" %in% object$par)  
    out$mu <- predict(object, what = "mu", type = type, terms = terms, se.fit = se.fit, ... )
    if ("sigma" %in% object$par)  
    out$sigma <- predict(object, what = "sigma", type = type, terms = terms, se.fit = se.fit, ... )
    if (  "nu" %in% object$par)     
    out$nu <- predict(object, what = "nu", type = type, terms = terms, se.fit = se.fit, ... )
    if ( "tau" %in% object$par)    
    out$tau <- predict(object, what = "tau", type = type, terms = terms, se.fit = se.fit, ... )
    if (fname%in%gamlss.bi.list)  out$bd <- object$bd  
if (output=="list")# LIST
   {# this part did  not working with gamlssML() so I have chanded lpred() with predict()
    attr(out, "family") <- object$family
    return(out)
   } # end list
if (output=="data.frame")# data.frame
 {
  df.out = as.data.frame(out)
  attr(df.out, "family") <- object$family
  return(df.out)
 }
if (output=="matrix")# MATRIX 
 {
  m.out = as.matrix(as.data.frame(out))
  attr(m.out, "family") <- object$family
  return(m.out)
 } 
}# OLD DATA ENDS  
else #  --------   if new data ----------------------------------------- 
{ # NEW DATA  CASE (ii)
     pN <- dim(newdata)[1]
 if (!(inherits(newdata, "data.frame")))
      stop("newdata must be a data frame ") # or a frame number
## NEW DATA CASE (ii) (a) se.fit=FALSE and use.weights==FALSE)
if ((use.weights==FALSE)&&(se.fit==FALSE))#  
 {    # use idividual predict.gamlss2()
  # if (output=="list") ## list
  # {
            out <- list()
  whetherFitted <- as.gamlss.family(object$family[1])$par
# if <par>.fix exists then set whetherFitted FALSE ------------------------
  whetherFitted <- as.list(unlist(whetherFitted) &
                             !paste0(names(whetherFitted), ".fix") %in% names(object))  
    if ("mu" %in% object$par)
        out$mu  <- if (whetherFitted$mu) 
          predict(object,newdata=newdata, what = "mu", type = type, 
                  terms = terms, data=DatA ) 
                  else rep(fitted(object)[1], pN)
    if ("sigma" %in% object$par)
      out$sigma <- if (whetherFitted$sigma) 
           predict(object, newdata=newdata, data=DatA,   what = "sigma",
                   type = type, terms = terms) 
                   else  rep(fitted(object, "sigma")[1], pN)
    if ("nu" %in% object$par) 
         out$nu <- if (whetherFitted$nu) predict(object, newdata=newdata, 
                    data=DatA,  what = "nu", 
                    type = type, terms = terms) 
                   else  rep(fitted(object, "nu")[1], pN)
    if ("tau" %in% object$par)  
        out$tau <- if (whetherFitted$tau) predict(object, newdata=newdata, 
                    data=DatA , what = "tau", type = type, terms = terms)  
                    else  rep(fitted(object, "tau")[1], pN)
    #if (fname%in%gamlss.bi.list)
if (any(ResCha%in%names(newdata))) 
    {
      TheY <-  newdata[,ResCha]
      if (fname%in%gamlss.bi.list)
      {
        if(is.matrix(TheY))
        {
          out$y <- TheY[,1]
          out$bd <- rowSums(TheY) 
        } else
        {
          out$y <- TheY
          out$bd <- rep(1, length(TheY)) 
        }  
      } else 
      {
        out$y <- TheY
      }  
     } # if list  finish here 
  if   (output=="list")#  list
  {
    attr(out, "family") <- object$family
    return(out)  
  }
  if (output=="data.frame")# data.frame
   {
     df.out = as.data.frame(out)
     attr(df.out, "family") <- object$family
     return(df.out)
   }
   if (output=="matrix")# MATRIX 
   {
     m.out = as.matrix(as.data.frame(out))
     attr(m.out, "family") <- object$family
     return(m.out)
   }
#######################################################################   
} # end  CASE (ii) a  
#######################################################################     
# CASE (ii) (b) if use.weights==TRUE) or se.fit==TRUE
####################################################################### 
if ((use.weights==TRUE)||(se.fit==TRUE))
 {
## here we should check whether newdata has ResCha
         ifY <- any(ResCha%in%names(newdata))     
##       keep only the same variables 
##       this assumes that all the relevant variables will be in newdata
##       what happens if not?     
##      The question is if newdata does not contain all variables what we do 
        DatA <- DatA[match(names(newdata),names(DatA), nomatch = 0)]# mikis 
## merge the two data together
        DatA <- concat(DatA,newdata)
## get the old weights
  oldweights <- object$weights
## create the new weights    
          pN <- dim(newdata)[1]
DatA$Weights <- c( rep(1, N)*oldweights, rep(set.to, pN)) 
## get the old response 
if  (!ifY) 
{
  oldresponse <- object$y
 yval <-  if (!is.numeric(y.value)) rep(do.call(y.value, list(x = oldresponse)), pN)
          else   rep(y.value, pN)
 DatA[[ResCha]] <- c(oldresponse, yval)
}
## now ifY is true then yvar should take values from the newdata
          # else{ # otherwise use y.value
          #   if (!is.numeric(y.value))  rep(do.call(y.value, list(x=oldresponse)), pN)
          #   else rep(y.value, pN)  
          # } 
##   cat("y", yval)
##     DatA[[ResCha ]] <- c(oldresponse, yval)
## refitting the model
     newobj <-  update(object, data=DatA, weights=DatA$Weights, trace=FALSE)
## needs some checking here   
     DiffAIC  <- AIC(object, newobj, k=0)
if (abs(diff(DiffAIC[,2]))>1)   warning("the global deviances differ more that 1 ")
if (abs(diff(DiffAIC[,1]))>.5)  warning("the degrees of freedom differ more that 0.5 ")
cat("dev:",deviance(object),deviance(newobj), deviance(object)-deviance(newobj),"\n")
cat("df:",object$df.fit,newobj$df.fit, object$df.fit-newobj$df.fit,"\n")
## prepare the output
if (fname%in%gamlss.bi.list) # if binomial
{
            prematrix <- matrix(0, ncol=length(object$parameters)+2, nrow=pN)
  colnames(prematrix) <- c("y", object$parameters, "bd")
} else 
{
  prematrix <- matrix(0, ncol=length(object$parameters)+lRes, nrow=pN)
  colnames(prematrix) <- c(ResCha, object$parameters)
   prematrix[,ResCha] <- as.matrix(newdata[,ResCha ])
}  
# binomial
if (fname%in%gamlss.bi.list)  prematrix[, "bd"] <-  newobj$bd[DatA$source!="DatA"]   
if ( se.fit)
     { 
               predse <- matrix(0, ncol=length(object$parameters), nrow=pN) 
     colnames(predse) <- object$parameters
for (i in object$parameters)
       { 
                PVaSE <- lpred(newobj, parameter=i, type=type, se.fit=TRUE)
        prematrix[,i] <- PVaSE$fit[(N+1):(N+pN)]  
           predse[,i] <- PVaSE$se.fit[(N+1):(N+pN)]  
       }
     colnames(predse) <- paste(colnames(predse), ".se", sep="")
              listout <- if (output=="matrix") cbind(prematrix, predse) else as.list(data.frame(cbind(prematrix, predse)))
      }# else
#      { 
# for (i in object$parameters)
#      { prematrix[,i] <- lpred(newobj, parameter=i, type=type)[(N+1):(N+pN)]  }
#              listout <- if (output=="matrix") prematrix else as.list(data.frame(prematrix))  
#        #  lapply(seq_len(ncol(prematrix)), function(i) prematrix[,i]) 
#      }    
}## EBD CASE (ii) b 
  if   (output=="list")#  list
     {
       attr(listout, "family") <- object$family
       attr(listout, "response") <- ResCha
       return(listout)  
     }
     if (output=="data.frame")# data.frame
     {
       df.out = as.data.frame(listout)
       attr(df.out, "family") <- object$family
       attr(df.out, "response") <- ResCha
       return(df.out)
     }
     if (output=="matrix")# MATRIX 
     {
       m.out = as.matrix(as.data.frame(listout))
       attr(m.out, "family") <- object$family
       attr(m.out, "response") <- ResCha
       return(m.out)
     }
} ## END CASE (ii)
} ## end of the function
################################################################################
################################################################################
################################################################################
################################################################################
#-------------------------------------------------------------------------------
################################################################################
################################################################################
################################################################################
################################################################################
