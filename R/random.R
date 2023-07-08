#-------------------------------------------------------------------------------
# this function is a modifyied version of  Hastie's S-plus gam test function
# it allows a random effect fit for a factor
# last modified August 2019
# TO DO:  i) check for weighted out obseervations
#        ii) put prediction OK 
#       iii) check if weighted levels 
random <- function(x, df = NULL, lambda = NULL, start=10) 
{
  scall <- deparse(sys.call(), width.cutoff = 500L)
 if(!inherits(x, "factor")) # | !is.category(xvar))
    stop("random() expects a factor as its first argument")
 if (!is.null(df))
   {
     nlevelsM <- nlevels(x)
     df <- if(df>=nlevelsM) nlevelsM  else df
     df <- if(df<=0) 0.001  else df
 }
   xvar <-  xvar <- rep(0,  length(x)) # C(x, rep(0, length(levels(x))), 1) # puts zero in the X matrix 
 attr(xvar, "call") <- substitute(gamlss.random(data[[scall]], z, w, ...))
 attr(xvar, "df") <- df
 attr(xvar, "lambda") <- lambda
 attr(xvar, "start") <- start
 attr(xvar, "contrasts") <- C(x, rep(0, length(levels(x))), 1) 
 class(xvar) <- c("smooth", class(xvar))
 xvar
}
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# this function is take from Hastie's gam 
# last change MS Thursday, April 25, 2002 at 13:31
gamlss.random <- function(x, y, w, xeval = NULL, ...) 
{
# local function----------------------------------------------------------------  
    df.inv <- function(n, df, lambda = sum(n)/df - mean(n), iterations = 10)
       { # given df find lambda
        if(df > length(n))
            return(0)
        current.df <- sum(n/(n + lambda))
        if(abs((df - current.df)/df) < 0.0001 | iterations == 1)
            lambda
        else {
            lambda <- exp(log(lambda) + (current.df - df)/(sum((n * lambda)/(n +lambda)^2)))
            Recall(n, df, lambda, iterations - 1)
             }
        }
#-------------------------------------------------------------------------------
regpen <- function(y, w, x, nw, lambda)
{
  beta <- tapply(w * y, x, sum)/(nw + lambda)
    df <- sum(nw[non.zero]/(nw[non.zero] + lambda))
   out <- list(beta=beta, edf=df)     
}
#-------------------------------------------------------------------------------
if (is.null(xeval))
{
          df <- as.vector(attr(x,"df"))
      lambda <- as.vector(attr(x, "lambda"))
       start <- as.vector(attr(x, "start"))
          xx <-  as.vector(attr(x, "contrasts"))
          nw <- tapply(w, xx, sum) 
    non.zero <- !is.na(nw)
        sig2 <- tau2 <- 0
           N <- sum(w!=0) # DS+FDB 3-2-14   
           n <- length(xx)
    # if lamnda is not 
# case 1 if lambda is know    
if (is.null(df)&&!is.null(lambda)||!is.null(df)&&!is.null(lambda))
 {
          fit <- regpen(y, w, xx, nw, lambda)
           fv <- fit$beta[xx]        
 } 
if (is.null(df)&&is.null(lambda)) # case 2 if both lambda and df are NULL
 {
   lambda <- start
   for (it in 1:50) 
         {
           fit  <- regpen(y, w, xx, nw, lambda) # fit model
         gamma. <- fit$beta  # get the gamma differences
             fv <-  fit$beta[xx]            # fitted values
           sig2 <- sum(w * (y - fv) ^ 2) / (N - fit$edf)
           tau2 <- sum(gamma. ^ 2) / (fit$edf)# Monday, March 16, 2009 at 20:00 see LNP page 279
     lambda.old <- lambda
         lambda <- sig2 / tau2 # maybe only 1/tau2 will do since it gives exactly the EM results see LM-1
     if (lambda<1.0e-7) lambda<-1.0e-7 # DS Saturday, April 11, 2009 at 14:18
         # cat("iter tau2 sig2",it,tau2, sig2, '\n')
     if (abs(lambda-lambda.old) < 1.0e-7||lambda>1.0e7) break
     # assign(startLambdaName, lambda, envir=gamlss.env)
     #cat("lambda",lambda, '\n')
         }
 }
if (!is.null(df)&&is.null(lambda)) # case 3 : if df are required----------------
 {
  lambda <- df.inv(nw[non.zero], df)
   fit <- regpen(y, w, xx, nw, lambda)
 }
  #  if(is.null(df))     df <- sum(non.zero)
  #  if(lambda == 0) lambda <- df.inv(nw[non.zero], df)
  #   df <- sum(nw[non.zero]/(nw[non.zero] + lambda))
     beta <- fit$beta #tapply(w * y, x, sum)/(nw + lambda)
       fv <- fit$beta[xx]
      var <- as.vector(w/(nw[xx] + lambda))
residuals <- as.vector(y -fv )
  coefSmo <- list(  coef = fit$beta,
                  lambda = lambda, 
                     edf = fit$edf, 
                   sigb2 = tau2, 
                   sige2 = sig2,
                    sigb = if (is.null(tau2)) NA else sqrt(tau2),
                    sige = if (is.null(sig2)) NA else sqrt(sig2),
                      fv = fv,  
                  factor = xx,
                      se = sqrt(var))
    class(coefSmo) <- "random"
  list(x = seq(along = nw), y = fv, residuals = residuals, var = var, 
       nl.df = fit$edf, lambda=lambda, coefSmo=coefSmo) # MS 
  }
  else 
  {  
            xx <- as.vector(attr(x, "contrasts"))
           xxx <- xx[(length(y)+1):length(xx)]
      position <- 0
         rexpr <- regexpr("predict.gamlss",sys.calls())
for (i in 1:length(rexpr)){ 
      position <- i 
if (rexpr[i]==1) break}
cat("New way of prediction in random()  (starting from GAMLSS version 5.0-6)", "\n")    
gamlss.environment <- sys.frame(position)
         param <- get("what", envir=gamlss.environment)
        object <- get("object", envir=gamlss.environment)
            TT <- get("TT", envir=gamlss.environment) 
 smooth.labels <- get("smooth.labels", envir=gamlss.environment)
    # needs to chect whether levels exist
           obj <- getSmo(object, parameter= param, which=which(smooth.labels==TT))
        FACTOR <- obj$factor      
            FV <- obj$fv 
    pred <- FV[xxx]  
  }
}
#---------------------------------------------------------------
plot.random <- function(x,...)
{
  plot(levels(x$factor), x$coef, type="h", xlab="levels", ylab="coefficients")
  abline(h=0)
}
#---------------------------------------------------------------
coef.random <- function(object, ...)
{
  as.vector(object$coef)
}
#--------------------------------------------------------------
fitted.random<- function(object, ...)
{
  as.vector(object$fv)
}
#--------------------------------------------------------------
print.random  <- function (x, digits = max(3, getOption("digits") - 3), ...) 
{   
  cat("Random effects fit using the gamlss function random() \n")
  cat("Degrees of Freedom for the fit :", x$edf, "\n")
  cat("Random effect parameter sigma_b:", format(signif(x$sigb)), "\n")  
  cat("Smoothing parameter lambda     :", format(signif(x$lambda)), "\n") 
}
#-------------------------------------------------------------
#-------------------------------------------------------------
getMarginal <- function(object, 
                  random.term = NULL, 
                    parameter = "mu",
                       method = c("integrate", "qfunction", "random", "none" ) 
                       )
{
if (!is.gamlss(object))  stop(paste("This is not a gamlss object", "\n", ""))
  method <- match.arg(method) 
if (is.null(random.term)) stop(paste("The random term is compulsory", "\n", ""))
# check id random.term exist in call 
#if (any(regexec(random.term,deparse(object$call))==TRUE))  
#     stop("term", random.term, "is not included in the fit" ,"\n", "")
# get the linear predictor
     lp <- eval(parse(text=(paste0("object$", parameter,".lp", sep=""))))
# get the contribution of the random term  
     rt <- eval(parse(text=paste0("object$", parameter,".s[,\"", random.term, "\"]" )))
# the linear predictor without e contribution of the random term 
   etam <- lp-rt
#   I need sigma_b
 sigma_b  <- eval(parse(text=paste0("getSmo(", deparse(substitute(object)), ", parameter=\"",parameter,"\")$sigb") ))
# I need te invesre link
      fam <- as.gamlss.family(substr(object$family[1], 1L, 30L)) # first family
 inv.link <- eval(parse(text=(paste("fam$", parameter,".linkinv", sep="")))) 
      out <- switch(method,
       "integrate"= sapply(etam, function(i){#Method 1: Integrating out random effects
                    integ <- function(x, pred, sd) inv.link(pred + x)*dNO(x, mu = 0, sigma = sigma_b)
                    integrate(integ, lower = -Inf, upper = Inf, pred = i)$value}),
       "qfunction"= sapply(etam, function(i, fam=fam){#Method 2
                    u <- qnorm(seq(0.001, 0.999, 0.001))*sigma_b
                    predu <- function(pred, u) {mean(inv.link(pred + u))}
                    predu(pred = i, u = u)}),
          "random"= sapply(etam, function(i){ 
                    u <- rNO(10000, mu=0, sigma=sigma_b)
                    predu2 <- function(pred, u) {mean(inv.link(pred + u))}
                    predu2(pred = i, u = u)}),
            "none"= inv.link(etam)
                      )
    out
}


