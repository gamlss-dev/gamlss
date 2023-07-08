####################################################################################
####################################################################################
####################################################################################
# this function is needed in ALE but it should be possibly use to all partial effect
#  functions 
#  for any use with pe_ we need to be able to create the list with all formulea from 
#  parameters
getAll.formulae2data <- function(obj, data)
{
  if (!is(obj,"gamlss")) stop("it needs a gamlss object")
  param <- obj$parameters
  lparam <- length(param)
  formul <- list()
  for (i in 1:lparam) 
  {
    # browser()
    #                           obj[[paste0(param[i],".formula")]][[3]][1:2]
    formul[[param[i]]] <- as.formula(obj[[paste0(param[i],".formula")]])
  }
  DaTa <-  Formulae2data(formul, data=data) #
  DaTa
}
################################################################################
################################################################################
################################################################################
Formulae2data <- function(formula = list(), data=NULL, weights=NULL, subset=NULL, 
                          na.action, print = TRUE  )
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
      f1 <- ff <- paste(paste(form[[2]],form[[1]]), deparse(form[[3]]), "+")
      # now add the of he formulae    
      for (i in 2:lenList)
      {
        ff <- if (i==lenList) paste(ff, deparse(formula[[i]][[2]]))
        else paste(ff, deparse(formula[[i]][[2]]),"+")
        
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
  if (!is.null(data)&&!inherits(data,"data.frame")) warning("data is not a data frame class attributes will be lost")
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
  all.vars
}
################################################################################
################################################################################
################################################################################