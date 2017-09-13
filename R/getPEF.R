###### Calculates the partial effect and elasticity of a continuous explanatory 
####### variable X.
###### It takes a GAMLSS object and for the range of a continous variable x, (and
###### by fixing the rest of the explanatory variables at specific values), 
###### calculates the effect x has on the specific distribution parameter (or its predictor).
###### The infuence function is then approximated using the splinefun() in R and saved.
###### The saved function and its derivatives can be used later. 
###### The fixed values for the rest of the explanatory variables are  
###### by default the median  
###### for continuous variables and the most commonly occuring level for 
###### factors 
getPEF <- function(obj = NULL, # the gamlss object
                  term = NULL, # which term to get the derivative
                  data = NULL, # which data is needed here
              n.points = 100,  # number of points needed for evaluating the function
             parameter = c("mu", "sigma", "nu", "tau"), # which parameter
                  type = c("response", "link"),
#                 pred = TRUE, #
                   how = c("median", "last"),
              fixed.at = list(), # see below (1)
                  plot = FALSE) # whether to plot
{
#  a named list of the values to use for the other predictor terms. Variables omitted from this list will have values set to the median for continuous variables and the most commonly occuring level for factors or the last observation
  if (is.null(obj)||!class(obj)[1]=="gamlss") stop("Supply a standard GAMLSS model in obj")
  if (is.null(term))  stop("The model term is not set")
            how <- match.arg(how)
           type <- match.arg(type)
      parameter <- match.arg(parameter)
if (any(grepl("data", names(obj$call)))) 
      {
        DaTa <- if (startsWith(as.character(obj$call["data"]), "na.omit"))
          eval(parse(text=as.character(obj$call["data"]))) else 
            get(as.character(obj$call["data"]))	
      }
      else if (is.null(data)) stop("The data argument is needed in obj")
            mat <- matrix(0, nrow = dim(DaTa)[1]+n.points, ncol =dim(DaTa)[2])
       dat.temp <- as.data.frame(mat)
names(dat.temp) <- v.names <- names(DaTa)
  #dat.temp$Date=obj$data$Date
            pos <- which(names(dat.temp)==term)
  if (pos<1) stop("supply a continuous term")
  if (is.factor(DaTa[,pos])) stop("the getPEF() is not suitable for factors")
## the values for which the function will be evaluated 
## we stick with min and max of x because we try to avoid trables with extrapolation            
           xvar <-  seq(from = min(DaTa[,pos]), to=max(DaTa[,pos]), length.out=n.points)
## creating new data         
  for (i in 1:dim(dat.temp)[2])
  {
    if(pos==i)                      # if the variable of interest
    {                               # new data for x is
   dat.temp[,i] <- c(DaTa[,i],xvar) # min(x) to max(x)
    }
    else                            # for all other variables
    {                               # if fixed.at is set gets priority
             ma <- fixed.at[[v.names[i]]]
      if (is.null(ma))              # if fixed.at in not set
      {
        if (how=="median")          # get the median for continuous 
                                    # or the level with high values for factor
        {
          ma <- if(is.factor(DaTa[,i])) levels(DaTa[,i])[which.max(table(DaTa[,i]))]
                else median(DaTa[,i])
        }
        if (how=="last")       # otherwise get the last values
        {
          ma <- if(is.factor(DaTa[,i])) levels(DaTa[,i])[which.max(table(DaTa[,i]))]
                else tail(DaTa[,i],1)
        }
      }
      dat.temp[,i] <- c(DaTa[,i],rep(ma,n.points))
    }
  } # end going thought the variables
## predict           
  fittted.orig <- predict(obj,newdata=tail(dat.temp,n.points), type = type, parameter = parameter)
  theFun  <- splinefun(xvar, fittted.orig)
  if (plot)
  {
   # browser() #def.par <- par(no.readonly = TRUE)
    layout(matrix(c(1,1,2,2), 2, 2, byrow = TRUE))
    plot(theFun(xvar)~xvar, ylab="s()", xlab=term, type="l")
    plot(theFun(xvar, deriv=1)~xvar, xlab=term,  ylab="ds/dx", type="l")
    abline(h=0)
    layout(1)
  }
  # theDeriv <-  splinefun(xvar, deriv)
     #deriv <- theFun(DaTa[,pos], deriv=1)
   #    out <- list(fittedValues = fittted.orig, theFun=theFun)
  invisible(theFun)
}
