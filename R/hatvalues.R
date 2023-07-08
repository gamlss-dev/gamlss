################################################################################
################################################################################
# the hat values for GAMLSS
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
# the linear leverage of a gamlss object
hatvalues.gamlss <- function(model, ...)
{
################################################################################
# local
################################################################################  
# Formulae2one <- function(formula, sigma=~1, nu=~1, tau=~1, data )
#   {
#     form <- formula(formula)
#     nform <- paste(paste(form[[2]],form[[1]]), deparse(form[[3]]), "+",  
#                    deparse(sigma[[2]]),"+",
#                    deparse(nu[[2]]),"+",
#                    deparse(tau[[2]]))[1]
#     ff<- formula(paste(nform, collapse = " "))
#     environment(ff) <- globalenv()
#     ff
# }
#####################################################################3########## 
  formulae2data <- function(formula = list(), data=NULL, weights=NULL, 
                            subset=NULL, na.action, print = TRUE  )
  {
    if (is(formula,"list"))
    {
      lenList <- length(formula)
      if (lenList==0) stop("no formula detected")
      if (lenList==1) 
      {
        ff <- deparse(formula[[1]])
      } else
      {
        # the first formula  
        form <- formula(formula[[1]])
        # create y~x+   
        f1 <- ff <- paste(paste(form[[2]],form[[1]]), 
                          deparse(form[[3]], width.cutoff = 500L), "+")
        # now add the of he formulae    
        for (i in 2:lenList)
        {
    ff <- if (i==lenList) paste(ff, deparse(formula[[i]][[2]], 
                                            width.cutoff = 500L))
          else paste(ff, deparse(formula[[i]][[2]], width.cutoff = 500L),"+")
        } 
      }
    } else if (is(formula,"formula")) {ff  <- deparse(substitute(formula))}
    else stop("The formula argument should be a formula or a list") 
    if (!is.null(weights)) 
    {
      # formula(paste(ff[[3]], collapse = " "))
      ff <- paste(ff, deparse(substitute(weights)), sep="+")
      # ff[[3]] <- paste(ff[[3]],deparse(substitute(weights)), sep="+")
    }
 environment(ff) <- globalenv()    # do I need this
        all.vars <- get_all_vars(ff, data=data)
if (!is.null(data)&&!inherits(data,"data.frame")) 
      warning("data is not a data frame class attributes will be lost")
    M <- dim(all.vars)[1]
    ## subsetting             
    if (!is.null(subset)) {
      r <- if (!is.null(data))  eval(substitute(subset), data,  parent.frame())
      else eval(substitute(subset),  parent.frame())
      if (!is.logical(r)) stop("'subset' must be logical")
      all.vars <- all.vars[r,]
      M <- dim(all.vars)[1]
      if (print) cat( M, "observations left after subsetting \n" )           
    }
    # it need a futher warning here      N <- dim(all.vars)[1]  
    # na.omit   
    all.vars <- na.omit(all.vars)                             # clear NA's
    N <- dim(all.vars)[1]     
    if (print) {if (M-N > 0) cat(M-N, "rows with NAs are deleted", "\n" )}
    if (print) cat( N, "observations with", dim(all.vars)[2], "variables \n")    
    attr(all.vars, "formula") <- ff
    return(all.vars)
  }
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
# getAll.formulae2data <- function(model, data)
# {
#   if (!is(model,"gamlss")) stop("it needs a gamlss object")
#   param <- model$parameters
#   lparam <- length(param)
#   formul <- list()
#   for (i in 1:lparam) 
#   {
#     formul[[param[i]]] <- as.formula(model[[paste0(param[i],".formula")]])
#   }
#   DaTa <-  Formulae2data(formul, data=data) #
#   DaTa
# }
################################################################################
# This in combination with  formulae2data()
# takes a gamlss object, extract its formulae and then 
# creates a data frame 
model2formulae2data <- function(model, data)
{
  if (!is(model,"gamlss")) stop("it needs a gamlss object")
  param <- model$parameters
  lparam <- length(param)
  formul <- list()
  for (i in 1:lparam) 
  {
    formul[[param[i]]] <- as.formula(model[[paste0(param[i],".formula")]])
  }
  DaTa <-  formulae2data(formul, data=data) #
  DaTa
}
################################################################################
################################################################################
# hatvalue.gamlss function starts here 
if (any(grepl("data", names(model$call)))) 
{
  DaTa <- if (startsWith(as.character(model$call["data"]), "na.omit")) 
             eval(parse(text = as.character(model$call["data"])))
         else get(as.character(model$call["data"]))
}
else if (is.null(data)) 
  stop("The data argument is needed in obj")
     weights <- model$weights
reduced_data <- model2formulae2data(model, data=DaTa)
        form <- attr(reduced_data, "formula")
reduced_data <- cbind(reduced_data, weights )
lev <- hatvalues(lm(form, weights=weights, data=reduced_data))
lev
}
################################################################################
################################################################################
################################################################################
################################################################################  
