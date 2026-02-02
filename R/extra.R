##---------------------------------------------------------------------------------------
## those are the extra functions needed for gamlss
## print.gamlss()
## family() as.family() family.default()
## gamlss.family()  gamlss.family.default()  as.gamlss.family()
## refit() refit.gamlss()  
## fitted.gamlss() coef.gamlss() residuals.gamlss() deviance.gamlss() predict.gamlss()
## formula.gamlss() is.gamlss() 
## lp(), model.frame.gamlss()
## IC(), AIC(), GAIC() : general information criterion
## weighted.mean() !out  
## checklink()
################################################################################
#                          refit                                               
################################################################################
#refit <- function (object,...)# MS Friday, June 6, 2003 at 11:43
#UseMethod("refit")

#refit.gamlss <- function (object, ...)
refit <- function (object, ...)
{  
   if (!is.gamlss(object))  stop(paste("This is not a gamlss object", "\n", ""))
   call <- object$call
   if (is.null(call)) 
        stop("need an object with call component")
      call <- object$call
    extras <- match.call(expand.dots = FALSE)$...
   new.cycle <- 2*object$control$n.cyc
 #  itercontrol <- list(iter=object$iter, n.cyc=new.cycle)# taken out Monday, March 17, 2008
   extras <- list(start.from= match.call()$object, iter=object$iter, n.cyc=new.cycle)
    existing <- !is.na(match(names(extras), names(call)))
   for (a in names(extras)[existing]) call[[a]] <- extras[[a]]
        if (any(!existing)) 
           {
            call <- c(as.list(call), extras[!existing])
            call <- as.call(call)
           }
   a <- eval(call, parent.frame())
  # rm(itercontrol,envir=sys.frame(sys.parent()))# taken out Monday, March 17, 2008
   a
}
################################################################################
#                         fitted.gamlss
################################################################################
fitted.gamlss<-function (object, what = c("mu", "sigma", "nu", "tau"), parameter= NULL, ... ) 
{
  what <- if (!is.null(parameter))  {
    match.arg(parameter, choices=c("mu", "sigma", "nu", "tau"))} else  match.arg(what)
if (!what%in%object$par) stop(paste(what,"is not a parameter in the gamlss object","\n"))
x <-  if (is.null(object$na.action)) object[[paste(what,"fv",sep=".")]]
      else napredict(object$na.action, object[[paste(what,"fv",sep=".")]])
x
}
################################################################################
#                         coef.gamlss
################################################################################
coef.gamlss<-function (object, what = c("mu", "sigma", "nu", "tau"), parameter= NULL, ... ) 
{
  what <- if (!is.null(parameter))  {
    match.arg(parameter, choices=c("mu", "sigma", "nu", "tau"))} else  match.arg(what)
if (!what%in%object$par) stop(paste(what,"is not a parameter in the object","\n"))
x <- object[[paste(what,"coefficients",sep=".")]]
x
}
##------------------------------------------------------------
##============================================================
# this probably shoulf go in the main gamlss 
coefAll <- function(obj, deviance=FALSE, ...)
{
  out <- list()
  if ("mu" %in% obj$par) #
    out$mu <- coef(obj, "mu")
  if ("sigma" %in% obj$par)  
    out$sigma <- coef(obj, "sigma")
  if ("nu" %in% obj$par)  
    out$nu <- coef(obj, "nu")
  if ("tau" %in% obj$par)  
    out$tau <-  coef(obj, "tau")
  if (deviance) out$deviance <- deviance(obj)
  return(out)
}
##============================================================
##------------------------------------------------------------
##============================================================

################################################################################
#                 residual or resid.gamlss
################################################################################
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
residuals.gamlss<-function (object, what = c("z-scores", "mu", "sigma", "nu", "tau"), 
                            type=c("simple","weighted","partial"), terms = NULL, ...) 
{
## Possible residuals
##  I)  z-scores  i) simpe or weighted
##                       simple just take the original:  object$residuals
##               ii) weighted   
##       
##                   a) "zero ones"  or "frequencies" 
##                               a1) Continuous rep(object$residuals, w)
##                               a2) discrete
##                   b) other    warning + object$residuals
##                                                                    
##  II)           "simple",  "weighted",  "partial"                   
## 
##  if what mu sigma nu tau  (the question here is whether w is needed) 
## type   simple       
##       weighted
##       partial
type <- match.arg(type)
what <- match.arg(what)
   w <- object$weights
if(what=="z-scores")                            #      if z-scores  (I) 
 {
  if (all(w==1)) x <- object$residuals          #      i) simpe or weighted
                                                #     ii) weighted   
  else if (all(trunc(w)==w))                    # a) "zero ones"  or  "frequencies" 
        { if (object$type== "Continuous")     x <- rep(object$residuals, w) # a1 continuous
          else{                                     # a2  discrete case
          y  <- rep(object$y, w)
      if ("mu"%in%object$parameters)        mu <- rep(fitted(object, "mu"),w)
      if ("sigma"%in%object$parameters)  sigma <- rep(fitted(object,"sigma"),w)
      if ("nu"%in%object$parameters)        nu <- rep(fitted(object,"nu"),w)
      if ("tau"%in%object$parameters)      tau <- rep(fitted(object,"tau"),w)
      if(any(object$family%in%.gamlss.bi.list)){ bd <- rep(object$bd,w)} # MS Wednesday, July 23, 2003 at 12:03   
          x <- eval(object$rqres)
              }  
         }  # now weights NOT "zero ones"  or  "frequencies"                             
else { warning("weights which are not frequencies are used: residuals remain unweighted")
           x <- object$residuals
     } # naresid(object$na.action, object$residuals)
 } # if not z-scores
 else   # II) If not z-scores
 { 
 if (!what%in%object$par) stop(paste(what,"is not a parameter in the object","\n"))
    wv  <- object[[paste(what, "wv", sep=".")]]
    l.p <- object[[paste(what, "lp", sep=".")]]
    wt  <- object[[paste(what, "wt", sep=".")]]
    os  <- object[[paste(what, "offset", sep=".")]]# offset is added 9-12-14 MS
    x   <- if (type == "simple") wv - l.p + os
         else if (type == "weighted") sqrt(wt)*(wv-l.p + os)
         else  (wv-l.p + os) +lpred(object, what = what, type = "terms", terms =terms) 
 } 
x
}
################################################################################

################################################################################
#                     deviance.gamlss
################################################################################
deviance.gamlss<-function(object, what = c("G", "P"), ...)
{ 
  what <- match.arg(what)
  if (what=="G")  x <- object$G.deviance
  else if (what=="P")  x <- object$P.deviance
  else stop("put G for Global or P for Penalized deviance")
 x 
}
################################################################################
#                   lp
################################################################################
## lm  see also lpred()
lp <-function (obj, what = c("mu", "sigma", "nu", "tau"), parameter= NULL, ... ) 
{
  what <- if (!is.null(parameter))  {
    match.arg(parameter, choices=c("mu", "sigma", "nu", "tau"))} else  match.arg(what)
if (!what%in%obj$par) stop(paste(what,"is not a parameter in the object","\n"))
x <- obj[[paste(what,"lp",sep=".")]]
x
}
################################################################################
#                   fv
################################################################################
## fv  see also fitted and lpred()
fv <-function (obj, what = c("mu", "sigma", "nu", "tau"), parameter= NULL, ... ) 
{
  what <- if (!is.null(parameter))  {
    match.arg(parameter, choices=c("mu", "sigma", "nu", "tau"))} else  match.arg(what)
if (!what%in%obj$par) stop(paste(what,"is not a parameter in the object","\n"))
x <- obj[[paste(what,"fv",sep=".")]]
x
}
################################################################################
#                 model.frame.gamlss
################################################################################
# new MS Thursday, June 24, 2004 at 13:45
model.frame.gamlss <-function (formula, what = c("mu", "sigma", "nu", "tau"), parameter= NULL, ...) 
{
  object <- formula
    dots <- list(...)
    what <- if (!is.null(parameter))  {
    match.arg(parameter, choices=c("mu", "sigma", "nu", "tau"))} else  match.arg(what)
    Call <- object$call
     parform <- formula(object, what)
    #parform <- object[[paste(what, "formula", sep=".")]]
    data <- if (!is.null(Call$data))  eval(Call$data, environment(formula$mu.terms))# MS 11-318 after suggestion in 
            else environment(formula$terms)
   Terms <- terms(parform)
      mf <- model.frame(Terms, data, xlev = object[[paste(what,"xlevels",sep=".")]])
   mf
}
################################################################################
#                 terms.gamlss
################################################################################
terms.gamlss <- function (x, what = c("mu", "sigma", "nu", "tau"), parameter= NULL, ...) 
{
  what <- if (!is.null(parameter))  {
    match.arg(parameter, choices=c("mu", "sigma", "nu", "tau"))} else  match.arg(what)
if (!what%in%x$par) stop(paste(what,"is not a parameter in the object","\n"))
    v <- x[[paste(what,"terms",sep=".")]]
    if (is.null(v)) 
        stop("no terms component")
    return(v)
}
################################################################################
#                  model.matrix
################################################################################
model.matrix.gamlss <- function (object, what = c("mu", "sigma", "nu", "tau"), parameter= NULL, ...) 
{
  what <- if (!is.null(parameter))  {
    match.arg(parameter, choices=c("mu", "sigma", "nu", "tau"))} else  match.arg(what)
if (!what%in%object$par) stop(paste(what,"is not a parameter in the object","\n"))
    v <- object[[paste(what,"x",sep=".")]]
    if (is.null(v)) 
        stop("no terms component")
    return(v)
}
################################################################################
#                 formula.gamlss
################################################################################
# formula.gamlss<-function (x, what = c("mu", "sigma", "nu", "tau"), parameter= NULL, ... ) 
# {
#   what <- if (!is.null(parameter))  {
#     match.arg(parameter, choices=c("mu", "sigma", "nu", "tau"))} else  match.arg(what)
#  if (!what%in%x$par) stop(paste(what,"is not a parameter in the object","\n")) 
#     fo <- x[[paste(what,"formula",sep=".")]]
#  ## the problem is when "." is in the formula, if true get formula from terms
#  if (length(fo)==2 && "."%in%strsplit(as.character(fo),split="")[[2]])# no resp 
#     fo <- formula(x[[paste(what,"terms",sep=".")]])
#  if (length(fo)==3 &&  "."%in%strsplit(as.character(fo), split = "")[[3]])# "."%in%strsplit(as.character(fo), split = "")[[3]]
#     fo <- formula(x[[paste(what,"terms",sep=".")]])
#     fo
# }
formula.gamlss<-function (x, what = c("mu", "sigma", "nu", "tau"), parameter= NULL, ... ) 
{
  what <- if (!is.null(parameter))  {
    match.arg(parameter, choices=c("mu", "sigma", "nu", "tau"))} else  match.arg(what)
  if (!what%in%x$par) stop(paste(what,"is not a parameter in the object","\n"))
  fo <- x[[paste(what,"formula",sep=".")]]
  # Making it not a "formula" to prevent incompatibility with package:formula.tools when using as.character()
  class(fo) = "formula_"
  ## the problem is when "." is in the formula, if true get formula from terms
  if ((length(fo)==2 && "."%in%strsplit(as.character(fo),split="")[[2]]) || (length(fo)==3 &&  "."%in%strsplit(as.character(fo), split = "")[[3]]))
    fo <- formula(x[[paste(what,"terms",sep=".")]])
  class(fo) = "formula"
  fo
}
################################################################################
#                 is.gamlss 
################################################################################
is.gamlss <- function (x) 
inherits(x, "gamlss")
#--------------------------------------------------------------------------------------
#is.formula <- function (x)
#inherits(x,"formula")
################################################################################
#                   IC
################################################################################
IC <- function(object, k=2)
{
if (is.gamlss(object)) object$G.dev+object$df.fit*k 
else stop(paste("this is not a gamlss object"))
}
################################################################################
#                 AIC.gamlss 
################################################################################
AIC.gamlss <- function (object, ..., k = 2, c = FALSE) 
{
 if (length(list(...))) 
      {
      object <- list(object, ...)
    
    isgamlss <- unlist(lapply(object, is.gamlss))
    if (!any(isgamlss)) stop("some of the objects are not gamlss")
          df <- as.numeric(lapply(object, function(x) x$df.fit))
           N <- as.numeric(lapply(object, function(x) x$N))
         Cor <- if ((k == 2)&&(c==TRUE)) (2*df*(df+1))/(N-df-1) else rep(0, length(object)) 
         AIC <- as.numeric(lapply(object, function(x) x$G.dev+x$df.fit*k ))+Cor  
         val <- as.data.frame(cbind(df,AIC))
        Call <- match.call()
      Call$k <- NULL
      Call$c <- NULL 
 row.names(val) <- as.character(Call[-1])
        val  <- val[order(AIC),]
        val
     }
  else 
     { val <- if (is.gamlss(object)) object$G.dev+object$df.fit*k 
                 else stop(paste("this is not a gamlss object"))    
              if ((k == 2)&&(c==TRUE)) val <- val + (2*object$df.fit*(object$df.fit+1))/(object$N-object$df.fit-1) 
       val 
      }
}
################################################################################
GAIC <- function(object, ..., k=2, c=FALSE)
{
  UseMethod("GAIC")
}
################################################################################
GAIC.default <- function(object, ...) {
  stop("No method for objects of class ", class(object))
}
################################################################################
GAIC.gamlss <- function(object,..., k = 2, c = FALSE ) #UseMethod("AIC")
{
 #GAIC<- function(object,..., k = 2, c = FALSE )
 if (length(list(...))) 
      {
      object <- list(object, ...)
    isgamlss <- unlist(lapply(object, is.gamlss))
    if (!any(isgamlss)) stop("some of the objects are not gamlss")
          df <- as.numeric(lapply(object, function(x) x$df.fit))
           N <- as.numeric(lapply(object, function(x) x$N))
         Cor <- if ((k == 2)&&(c==TRUE)) (2*df*(df+1))/(N-df-1) else rep(0, length(object)) 
         AIC <- as.numeric(lapply(object, function(x) x$G.dev+x$df.fit*k ))+Cor  
         val <- as.data.frame(cbind(df,AIC))
        Call <- match.call()
      Call$k <- NULL
      Call$c <- NULL 
 row.names(val) <- as.character(Call[-1])
        val  <-  val[order(AIC),]
        val
     }
  else 
     { val <- if (is.gamlss(object)) object$G.dev+object$df.fit*k 
                 else stop(paste("this is not a gamlss object"))
       if ((k == 2)&&(c==TRUE)) val <- val + (2*object$df.fit*(object$df.fit+1))/(object$N-object$df.fit-1) 
       val 
      }
}
################################################################################
################################################################################
# GAIC.table <- function(object,..., 
#                    k = c(2, 3.84, round(log(length(object$y)),2)),
#                    text.to.show=NULL) 
# {
#   .Defunct(GAIC_table, package="gamlss")
# }
################################################################################
GAIC_table <- function(object,..., 
                       k = c(2, 3.84, round(log(length(object$y)),2)),
                       text.to.show=NULL) 
{
  objects <- list(object, ...)
    norow <- length(objects)
       lk <- length(k)
      MM  <- matrix(0, nrow=norow, ncol=lk+1)
 isgamlss <- unlist(lapply(objects, is.gamlss))
  if (!any(isgamlss)) stop("some of the objects are not gamlss")
   MM[,1] <- unlist(lapply(objects, function(x) x$df.fit))
        N <- as.numeric(lapply(objects, function(x) x$N))
  for (j in 1:lk)
  {
    MM[,1+j] <- unlist(lapply(objects, function(x) x$G.dev+x$df.fit*k[j] ))
  }
           val <- as.data.frame(MM)
          Call <- match.call()
        Call$k <- NULL
row.names(val) <- if (is.null(text.to.show)) as.character(Call[-1])
                  else text.to.show
    names(val) <- c("df", paste("k=",as.character(k), sep=""))
  for (i in 1:length(k)) cat("minimum GAIC(k=",k[i],") model:", 
                             row.names(val)[which.min(unclass(val[1+i][[1]]))], "\n")
  val
}
################################################################################
################################################################################
# Creating AIC weights
# the weights can be interpreted as posterior probabilities
# i.e prob(model/data)
# but more like as evidence on how the diffrent models 
# compared with each other
#  diff.dev should be a function of the sample size
################################################################################
################################################################################
GAIC_scaled <- function(object,..., 
                      k = 2,      # which column
                      c = FALSE,  # whether corrected GAIC
                   plot = TRUE,   #
               text.cex = 0.7,    # the model/distribution text size
                  which = 1,      # for families which GAIC
               diff.dev = 1000,
           text.to.show = NULL,
                    col = NULL,   # color for the bars
                  horiz = FALSE)  #  verical or horizontal
{
if (is.matrix(object))
  {
    fnas <- !is.na(object[,which])
      dm <- dim(object)
# object
    wmin <- which.min(object[,which])  #apply(object, 2, which.min)
    dAIC <- (object[,which]-object[wmin,which])
    dAIC <- ifelse( dAIC > diff.dev, NA, dAIC)
    #    which.max(dAIC)
    oAIC <- 1-dAIC/dAIC[which.max(dAIC)] #
    #   }
  if (plot)
    {    
      which.na <- !is.na(oAIC)   
     namesDist <- rownames(object)[which.na] 
        scaled <- oAIC[which.na]
      names(scaled) <- ""
          ind  <- c(1:sum(which.na))
     if (horiz)
      {
        BP <- barplot(scaled, xlim=c(-0.25, 1.1), horiz=TRUE,  col=col )
        points(scaled, BP, pch=20)
        abline(v=0, col="black")
        abline(v=1, col="gray")
        grid()
        box()
        text(-0.01, BP, namesDist, srt = 0, cex=text.cex, adj=1) 
      } 
      else 
      { 
        BP <- barplot(scaled,  ylim=c(-0.25, 1.1), col=col)
        points(BP, scaled, pch=20)
        abline(h=0, col="gray")
        abline(h=1, col="gray")
        grid()
        box()
        text(BP,-0.01, namesDist, srt = 90, cex=text.cex, adj=1)
      }
    }
  } else
if (length(list(...))) 
    {
      object <- list(object, ...)
    isgamlss <- unlist(lapply(object, is.gamlss))
if (!any(isgamlss)) stop("some of the objects are not gamlss")
          df <- as.numeric(lapply(object, function(x) x$df.fit))
           N <- as.numeric(lapply(object, function(x) x$N))
         Cor <- if ((k == 2)&&(c==TRUE)) (2*df*(df+1))/(N-df-1) else rep(0, length(object)) 
         AIC <- as.numeric(lapply(object, function(x) x$G.dev+x$df.fit*k ))+Cor  
        dAIC <- (AIC-AIC[which.min(AIC)]) 
        oAIC <- 1-dAIC/dAIC[which.max(dAIC)]
         val <- as.data.frame(cbind(df,AIC, delta=round(dAIC,5), scaled=round(oAIC,4)))
        Call <- match.call()
      Call$k <- Call$c <- Call$plot <- Call$text.cex <- Call$which <- Call$diff.dev <- Call$horiz <- NULL 
      row.names(val) <- if (is.null(text.to.show)) as.character(Call[-1])
      else text.to.show
   if  (plot) 
    {
      ind <- 1:length(AIC)
      if (horiz)
        {
          BP <- barplot(oAIC,   xlim=c(-0.25, 1.1), horiz=TRUE,  col=col )
          points(oAIC,BP, pch=20)
          abline(v=0, col="black")
          abline(v=1, col="gray")
          grid()
          box()
          text(-0.01, BP, rownames(val), srt = 0, cex=text.cex, adj=1) 
        } else
        { 
          BP <- barplot(oAIC,   ylim=c(-0.25, 1.1),  col=col)
          points(BP, oAIC, pch=20)
          abline(h=0, col="gray")
          abline(h=1, col="gray")
          grid()
          box()
          text(BP,-0.01, rownames(val), srt = 90, cex=text.cex, adj=1)
        }
      }
      val  <-  val[,]
      val
    }
  else 
  { val <- if (is.gamlss(object)) object$G.dev+object$df.fit*k 
  else stop(paste("this is not a gamlss object"))
  if ((k == 2)&&(c==TRUE)) val <- val + (2*object$df.fit*(object$df.fit+1))/(object$N-object$df.fit-1) 
  val 
  }
}
################################################################################
################################################################################
# GAIC.weights <- function(object, ..., 
#                          k = 2,      # whic gAIC
#                          c = FALSE,  # whether corrected GAIC
#                       plot = TRUE,   #
#                   text.cex = 0.7,    # the model/distribution text size
#                      which = 1,      # for families which GAIC
#                   diff.dev = 1000)   # if difference i GD more that 1000 do not print
# {
# if (is.matrix(object))
#   {
#     fnas <- !is.na(object)
#     dm <- dim(object)
#     # object
#     wmin <- apply(object, 2, which.min)
#     dAIC <- wAIC <- matrix(0,nrow=dm[1], ncol=dm[2] )
#     for (i in 1:dm[2]){
#       dAIC[,i] <- (object[,i]-object[wmin[i],i])[fnas[i]]
#       dAIC[,i] <- ifelse( dAIC[,i]>diff.dev, NA, dAIC[,i])
#           dAIC <- 2*scale(dAIC) # this is rather abritary 
#       wAIC[,i] <- exp(-0.5*dAIC[,i])/sum(exp(-0.5*dAIC[,i]), na.rm=TRUE)
#     }
#     if (plot)
#     {    
#       which.na <- !is.na(wAIC[,which])   
#      namesDist <- rownames(object)[which.na] 
#        Weights <- wAIC[which.na,which]
#           ind  <- c(1:length(namesDist))
#       plot(Weights~ind, type="h", ylim=c(-0.22, max(Weights)+.05), xlab="families")
#       points(Weights, pch=20)
#       abline(h=0)
#       text(ind,-0.01, namesDist, srt = 90, cex=text.cex, adj=1) #
#     }
#   } else if (length(list(...))) 
#   {
#     object <- list(object, ...)
#   isgamlss <- unlist(lapply(object, is.gamlss))
# if (!any(isgamlss)) stop("some of the objects are not gamlss")
#         df <- as.numeric(lapply(object, function(x) x$df.fit))
#          N <- as.numeric(lapply(object, function(x) x$N))
#        Cor <- if ((k == 2)&&(c==TRUE)) (2*df*(df+1))/(N-df-1) else rep(0, length(object)) 
#        AIC <- as.numeric(lapply(object, function(x) x$G.dev+x$df.fit*k ))+Cor  
#       dAIC <- (AIC-AIC[which.min(AIC)]) 
#       #dAIC <- ifelse( dAIC > diff.dev, NA, dAIC)
#       sAIC <- 2*scale(dAIC) # this is rather abritary 
#       wAIC <- exp(-0.5*sAIC)/sum(exp(-0.5*sAIC))
#        val <- as.data.frame(cbind(df, AIC, delta=round(dAIC,5), weights=round(wAIC,4)))
#       Call <- match.call()
#     Call$k <- Call$c <- NULL 
# row.names(val) <- as.character(Call[-1])
#  if  (plot) 
#       {
#         ind <- 1:length(AIC)
#         plot(ind, wAIC , ylab="weights", xlab="models", 
#              type="h", col="black", ylim=c(-0.1, max(wAIC)+0.1))
#         points(1:length(AIC), wAIC, pch=20)
#         abline(h=0)
#         text(ind,-0.01, rownames(val), srt = 90, cex=text.cex, adj=1) 
#       }
#       val
#     }
#   else 
#   { val <- if (is.gamlss(object)) object$G.dev+object$df.fit*k 
#   else stop(paste("this is not a gamlss object"))
#   if ((k == 2)&&(c==TRUE)) val <- val + (2*object$df.fit*(object$df.fit+1))/(object$N-object$df.fit-1) 
#   val 
#   }
# }
################################################################################
################################################################################
##############################################################################
##############################################################################
## a small utility function to get the hat matrix from a weighted regression
## used in cs() and other additive terms 
## MS Tuesday, June 22, 2004 at 21:26
.hat.WX<-function (w,x) 
{
    p <- length(x)
    X  <- if (!is.matrix(x)) matrix(cbind(1,x), ncol=2)
          else x 
    k <- length(w)
    p <- dim(X)[1]
    if (p != k) 
        stop("`w' and 'x' are not having the same length")
    Is <- sqrt(w)
    if (any(!is.finite(Is))) 
        warning("diagonal weights has non-finite entries")
    WX <- X 
    WX[] <- Is * X
    h<-hat(qr(WX))
    h
}
################################################################################
################################################################################
# the generalised R-squared funtion
# This should work now with Binomial
Rsq <- function(object, type = c("Cox Snell","Cragg Uhler","both"),...)
{
  UseMethod("Rsq")
}
################################################################################
################################################################################
Rsq.default <- function(object, ...) {
  stop("No method for objects of class ", class(object))
}
################################################################################
################################################################################
Rsq.gamlss <- function(object, type = c("Cox Snell","Cragg Uhler","both"),...)
{
  type <- match.arg(type)
  if (!is.gamlss(object)) stop("this is design for gamlss objects only")
  #  m0 <- update(object,  formula=~1, sigma.formula=~1, nu.formula=~1, tau.formula=~1, trace=F)
  Y <- if (object$family[1]%in%.gamlss.bi.list) cbind(object$y, object$bd-object$y) else object$y
  suppressWarnings(m0 <- gamlssML(Y, family=object$family[1]))
  rsq1 <- 1-exp((2/object$N)*(logLik(m0)[1]-logLik(object)[1]))
  rsq2 <- rsq1/(1-exp((2/object$N)*logLik(m0)[1]))
  if (type=="Cox Snell") return(rsq1)
  if (type=="Cragg Uhler") return(rsq2)
  if (type=="both") return(list(CoxSnell=rsq1, CraggUhler=rsq2)) 
}
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# the function is in a new file now
# this intoducing the shifted log and logit links 
# MS Sunday, February 20, 2005 
# this will replace the make.link used for GLM's and GAM's
# the add link functions are 
# (1) lofshifted needing one parameters the left shift
# (2)logitshift.5 needing two parameters the left and rirght shift
# make.link.gamlss<-function (link, par = 1) 
#{
#    if (is.character(link) && length(grep("^power", link) > 0)) 
#        return(eval(parse(text = link)))
#    else if (!is.character(link) && !is.na(lambda <- as.numeric(link))) {
#        linkfun <- function(mu) mu^lambda
#        linkinv <- function(eta) pmax(.Machine$double.eps, eta^(1/lambda))
#         mu.eta <- function(eta) pmax(.Machine$double.eps, (1/lambda) * 
#                     eta^(1/lambda - 1))
#       valideta <- function(eta) all(eta > 0)
#    }
#    else switch(link, logit = {
#        linkfun <- function(mu) log(mu/(1 - mu))
#        linkinv <- function(eta) 
#            {
#            thresh <- -log(.Machine$double.eps)
#               eta <- pmin(thresh, pmax(eta, -thresh))
#                      exp(eta)/(1 + exp(eta))
#            }
#        mu.eta <- function(eta) 
#            {
#            thresh <- -log(.Machine$double.eps)
#               res <- rep(.Machine$double.eps, length(eta))
#            res[abs(eta) < thresh] <- (exp(eta)/(1 + exp(eta))^2)[abs(eta) < 
#                thresh]
#            res
#            }
#      valideta <- function(eta) TRUE
#    }, probit = {
#        linkfun <- function(mu) qnorm(mu)
#        linkinv <- function(eta) 
#           {
#            thresh <- -qnorm(.Machine$double.eps)
#               eta <- pmin(thresh, pmax(eta, -thresh))
#            pnorm(eta)
#           }
#        mu.eta <- function(eta) pmax(dnorm(eta), .Machine$double.eps)
#      valideta <- function(eta) TRUE
#    }, cloglog = {
#       linkfun <- function(mu) log(-log(1 - mu))
#       linkinv <- function(eta) pmax(.Machine$double.eps, pmin(1 - 
#                .Machine$double.eps, -expm1(-exp(eta))))
#        mu.eta <- function(eta) 
#            {
#           eta <- pmin(eta, 700)
#            pmax(.Machine$double.eps, exp(eta) * exp(-exp(eta)))
#            }
#      valideta <- function(eta) TRUE
#    }, identity = {
#       linkfun <- function(mu) mu
#       linkinv <- function(eta) eta
#        mu.eta <- function(eta) rep(1, length(eta))
#      valideta <- function(eta) TRUE
#    }, log = {
#       linkfun <- function(mu) log(mu)
#       linkinv <- function(eta) pmax(.Machine$double.eps, exp(eta))
#        mu.eta <- function(eta) pmax(.Machine$double.eps, exp(eta))
#      valideta <- function(eta) TRUE
#    }, sqrt = {
#       linkfun <- function(mu) mu^0.5
#       linkinv <- function(eta) eta^2
#        mu.eta <- function(eta) 2 * eta
#       valideta <- function(eta) all(eta > 0)
#    }, "1/mu^2" = {
#       linkfun <- function(mu) 1/mu^2
#       linkinv <- function(eta) 1/eta^0.5
#        mu.eta <- function(eta) -1/(2 * eta^1.5)
#      valideta <- function(eta) all(eta > 0)
#     }, logshiftto1 = {    # MS Saturday, February 19, 2005   
#       linkfun <- function(mu, shift = par, delta = .00001 )           
#                         log(mu-shift+delta)
#       linkinv <- function(eta, shift = par) shift+pmax(.Machine$double.eps, exp(eta)) 
#        mu.eta <- function(eta, shift = par, delta  = .00001) pmax(.Machine$double.eps, exp(eta))
#      valideta <- function(eta) TRUE   
#      },logitshift.5 = { # MS Saturday, February 19, 2005   
#       linkfun <- function(mu, shift = par )           
#                         log((mu-shift[1])/(shift[2]-mu))
#       linkinv <- function(eta,  shift = par) 
#            {
#            thresh <- -log(.Machine$double.eps)
#               eta <- pmin(thresh, pmax(eta, -thresh))
#                      shift[2]-(shift[2]-shift[1])/(1 + exp(eta))
#            } 
#        mu.eta <- function(eta, shift = par ) 
#            {
#            thresh <- -log(.Machine$double.eps)
#               res <- rep(.Machine$double.eps, length(eta))
#            res[abs(eta) < thresh] <- ((shift[2]-shift[1])*exp(eta)/(1 + exp(eta))^2)[abs(eta) < thresh]
#            res
#            }
#      valideta <- function(eta) TRUE       
#    }, inverse = {
#       linkfun <- function(mu) 1/mu
#       linkinv <- function(eta) 1/eta
#        mu.eta <- function(eta) -1/(eta^2)
#      valideta <- function(eta) all(eta != 0)
#    }, stop(paste(link, "link not recognised")))
#    list(linkfun = linkfun, linkinv = linkinv, mu.eta = mu.eta, 
#        valideta = valideta)
#}
##------------------------------------------------------------------------------
# numerical derivatives
# taken from the Writing R Extensions p 48
# modified to have a varied delta 
numeric.deriv <- function(expr, theta, delta = NULL,  rho=sys.frame(sys.parent()))
{
  eps <- sqrt(.Machine$double.eps)  
  ans <- eval(substitute(expr), rho)
 grad <- matrix(,length(ans), length(theta), dimnames=list(NULL, theta))
for (i in seq(along=theta)) 
  {
  old <- get(theta[i], envir=rho)
delta <-  if (is.null(delta)) eps * min(1, abs(old)) else delta
 assign(theta[i], old+delta, envir=rho)
 ans1 <- eval(substitute(expr), rho)
 assign(theta[i], old, envir=rho)
 grad[,i] <- (ans1 - ans)/delta
  }
attr(ans, "gradient") <- grad
ans
}
#-------------------------------------------------------------------------------
