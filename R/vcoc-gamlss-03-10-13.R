#  Mikis Stasinopoulos 22-03-13
# This is the new vcov.gamlss() function 
# it was created in Jan 2013 but impemented in gamlss on the 22 March 2013  
# it is using the function gen.likelihood

#-------------------------------------------------------------------------------
vcov.gamlss <- function (object, 
                           type = c("vcov", "cor", "se", "coef", "all"),
                         robust = FALSE, 
                ...) 
{
      type <- match.arg(type)
  if (!is.gamlss(object)) 
     stop(paste("This is not an gamlss object", "\n", ""))
  coefBeta <- list()
  for (i in object$par) 
  {
    if (i == "mu") 
      {
      if (!is.null(unlist(attr(terms(formula(object), specials = .gamlss.sm.list), 
                               "specials")))) 
        warning("addive terms exists in the mu formula standard errors for the linear terms maybe are not appropriate")
    }
    else 
    {
      if (!is.null(unlist(attr(terms(formula(object, i), 
                                     specials = .gamlss.sm.list), "specials")))) 
        warning(paste("addive terms exists in the ", 
                      i, " formula standard errors for the linear terms maybe are not appropriate"))
    }
  #    parname <- paste(i, "start", sep = ".")
     coefBeta <- c(coefBeta, coef(object, i))
  }
   betaCoef <- unlist(coefBeta)      
   like.fun <- gen.likelihood(object)
       hess <- optimHess(betaCoef, like.fun)
     varCov <- try(solve(hess))
      if (any(class(beta)%in%"try-error"))
      {stop("the Hessian matrix is singular probably the model is overparametrised")}
         se <- sqrt(diag(varCov))
       corr <- cov2cor(varCov) # cov/(se %o% se)
   coefBeta <- unlist(coefBeta)
  #names(coefBeta) <- attr(a$se, "names")
    if (robust)
    {
      K <- get.K(object)
      varCov <- varCov%*%K%*%varCov
      se <- sqrt(diag(varCov))
      corr <- cov2cor(varCov)
    }
  switch(type, vcov = varCov, cor = corr, se = se, coef = coefBeta, 
         all = list(coef = coefBeta, se = se, vcov = varCov, 
                    cor = corr))
}