#-----------------------------------------------------------------------------------
getSmo <- function(object, what = c("mu", "sigma", "nu", "tau"), which=1)
{
  if (!is.gamlss(object)) stop("this is design for gamlss objects only")
  what <- match.arg(what)
  if (!what%in%object$par) stop(paste(what,"is not a parameter in the object","\n")) 
  evalq(paste("object$",what,"[", which, "]",sep=""))
  AllSmo <- object[[paste(what,".coefSmo",sep="")]]
  Smo <- if (which==0) AllSmo else AllSmo[[which]]
  Smo
}
