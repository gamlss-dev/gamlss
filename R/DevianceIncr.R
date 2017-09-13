# deviance increment function
#---------------------------------------------------------------------
# the DevianceIncr() function will calculate the global deviance increment
# i)  for a fitted gamlss model or
# ii) for test data if  the option newdata is used
# large values for deviance increment indicate bad fit
# and for new data bad prediction
devianceIncr <- function(obj, newdata=NULL)
{
  if (!is.gamlss(obj)) stop("This works for gamlss objects only")
         pdf <- obj$family[1]
          fn <- eval(parse(text=pdf))()$G.dev.inc
        npar <- eval(parse(text=pdf))()$nopar
  parameters <- names(eval(parse(text=pdf))()$parameters)
  if ("mu"%in%parameters) 
          mu <- predict(obj, what = "mu", newdata = newdata, type = "response")
  if ("sigma"%in%parameters)
       sigma <- predict(obj, what = "sigma", newdata = newdata, type = "response")
  if ("nu"%in%parameters) 
          nu <- predict(obj, what = "nu", newdata = newdata, type = "response", 
                       na.action = na.omit)
  if ("tau"%in%parameters)
         tau <- predict(obj, what = "tau", newdata = newdata, type = "response", 
                       na.action = na.omit)
  if (is.null(newdata))
  {
    DINC <- switch(npar, fn(obj$y, mu=mu),
                 fn(obj$y, mu=mu, sigma=sigma),
                 fn(obj$y, mu=mu, sigma=sigma, nu=nu),
                 fn(obj$y, mu=mu, sigma=sigma, nu=nu, tau=tau))
  } else
  {
      YY <- eval(attributes(terms(formula(obj)))$variables[[2]], envir=newdata)
    DINC <- switch(npar, fn(obj$y, mu=mu),
                   fn(YY, mu=mu, sigma=sigma),
                   fn(YY, mu=mu, sigma=sigma, nu=nu),
                   fn(YY, mu=mu, sigma=sigma, nu=nu, tau=tau)) 
  }  
  DINC
}
#-----------------------------------------------------------------

