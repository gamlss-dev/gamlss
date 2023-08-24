\name{fitDist}
\alias{fitDist}
\alias{fitDistPred}
\alias{chooseDist}
\alias{chooseDistPred}
\alias{getOrder}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Fitting Different Parametric \code{gamlss.family} Distributions.
}
\description{
The function \code{fitDist()} is using the function \code{gamlssML()} to fit all relevant parametric \code{gamlss.family} distributions, specified by the argument \code{type}), to a single data vector (with no explanatory variables).  The final marginal distribution  is the one selected by the generalised Akaike information criterion with penalty \code{k}. The default is  \code{k=2} i.e AIC.

The function \code{fitDistPred()} is using the function \code{gamlssMLpred()} to fit all relevant (marginal) parametric \code{gamlss.family} distributions to a single data vector (similar to \code{fitDist()})  but the final model is selected by the minimum prediction global deviance. The user has to specify the training  and validation/test samples.


The function \code{chooseDist()} is using the function \code{update.gamlss()} to fit all relevant parametric (conditional)  \code{gamlss.family} distributions to a given fitted \code{gamlss} model.  The output of the function is a matrix with  rows the different distributions (from the argument \code{type}) and columns the different GAIC's (). The default argument for \code{k} are 2, for AIC, 3.84, for Chi square, and log(n) for BIC.  No final model is given by the function  like for example  in \code{fitDist()}.   The function \code{getOrder()} can be used to rank the columns of the resulting table (matrix).
The final model can be refitted using \code{update()}, see the examples.  

}
\usage{
fitDist(y, k = 2, 
    type = c("realAll", "realline", "realplus", "real0to1", "counts", "binom"), 
    try.gamlss = FALSE, extra = NULL, data = NULL,trace = FALSE, ...)

fitDistPred(y, 
    type = c("realAll", "realline", "realplus", "real0to1", "counts", "binom"), 
    try.gamlss = FALSE, extra = NULL, data = NULL, rand = NULL,
    newdata = NULL, trace = FALSE, ...)    
      
chooseDist(object, k = c(2, 3.84, round(log(length(object$y)), 2)), type = 
    c("realAll", "realline", "realplus", "real0to1", "counts", "binom","extra"), 
    extra = NULL, trace = FALSE, 
    parallel = c("no", "multicore", "snow"), ncpus = 1L, cl = NULL, ...)

chooseDistPred(object, type = c("realAll", "realline", "realplus", 
     "real0to1", "counts", "binom", "extra"), extra = NULL, 
     trace = FALSE, parallel = c("no", "multicore", "snow"), 
     ncpus = 1L, cl = NULL, newdata = NULL, rand = NULL, ...)
      
getOrder(obj, column = 1)      
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{the data vector}
  \item{object, obj}{a GAMLSS fitted model}
  \item{k}{the penalty for the GAIC with default values \code{k=2} the standard AIC.  In the case of the function \code{chooseDist()} \code{k} can be a vector i.e. \code{k= c(2, 4, 6)} so more than one GAIC are saved.}
  \item{type}{the type of distribution to be tried see details}
  \item{try.gamlss}{this applies to functions \code{fitDist()} and \code{fitDistPred()}.  It allows if  \code{gamlssML()} fail to fit the model to try \code{gamlss} instead. This will slow up things for big data.}
  \item{extra}{whether extra distributions should be tried, which are not in the \code{type} list. Note that the function \code{chooseDist()} allows the fitting of only the `extra' distributions. This can be achieved if \code{extra} is set i.e. \code{extra=c("GA", "IG", "GG")} and type is set to extra i.e. \code{type="extra"}. }
  \item{data}{the data frame where \code{y} can be found, only for functions \code{fitDist()} and \code{fitDistPred()} }
   \item{rand}{For \code{fitDistPred()} a factor with values 1 (for fitting) and 2 (for predicting).}
  \item{newdata}{The prediction data  set (validation or test).}
  \item{trace}{whether to print during fitting. Note that when \code{parallel}  is 'multocore' or "snow" \code{"trace"} is not  produce any output.}
    \item{parallel}{The type of parallel operation to be used (if any). If missing, the default is "no".}
   \item{ncpus}{integer: number of processes to be used in parallel operation: typically     one would chose this to the number of available CPUs.}
  \item{cl}{This is useful for snow clusters, i.e. \code{parallel = "snow"}, when the clusters are created in advance. If not supplied, a cluster on the local machine is created for the duration of the call.} 
  \item{column}{which column of the output matrix to be ordered according to best GAIC}
  \item{\dots}{for extra arguments to be passed to gamlssML() to gamlss()}
}
\details{
The following are the different \code{type} argument:
\itemize{
  \item{\code{realAll}:}{ All the \code{gamlss.family} continuous distributions defined on the real line, i.e. \code{realline} and the real positive line i.e. \code{realplus}}
  
  \item{\code{realline}:}{ The \code{gamlss.family} continuous distributions : "NO", "GU", "RG" ,"LO", "NET", "TF", "TF2", "PE","PE2", "SN1", "SN2", "exGAUS", "SHASH", "SHASHo","SHASHo2", "EGB2", "JSU", "JSUo", "SEP1", "SEP2", "SEP3", "SEP4",  "ST1", "ST2", "ST3", "ST4", "ST5", "SST", "GT"}
  
  \item{\code{realplus}:}{ The \code{gamlss.family} continuous distributions in the positive real line: "EXP", "GA","IG","LOGNO", "LOGNO2","WEI", "WEI2", "WEI3", "IGAMMA","PARETO2", "PARETO2o", "GP", "BCCG", "BCCGo", "exGAUS", "GG", "GIG", "LNO","BCTo", "BCT", "BCPEo", "BCPE", "GB2"}
  
  \item{\code{real0to1}:}{ The \code{gamlss.family} continuous distributions from 0 to 1: "BE", "BEo", "BEINF0", "BEINF1", "BEOI", "BEZI", "BEINF", "GB1""}
  
 \item{\code{counts}:}{ The \code{gamlss.family} distributions for counts: "PO", "GEOM", "GEOMo","LG", "YULE", "ZIPF", "WARING", "GPO", "DPO", "BNB", "NBF","NBI", "NBII", "PIG", "ZIP","ZIP2", "ZAP", "ZALG", "DEL", "ZAZIPF", "SI", "SICHEL","ZANBI",  "ZAPIG", "ZINBI",  "ZIPIG", "ZINBF", "ZABNB", "ZASICHEL", "ZINBF",  "ZIBNB", 
 "ZISICHEL"}

\item{\code{binom}:}{ The \code{gamlss.family} distributions for binomial type data :"BI", "BB", "DB", "ZIBI", "ZIBB", "ZABI", "ZABB"}

The function \code{fitDist()} uses the function \code{gamlssML()} to fit the different models, the function  \code{fitDistPred()} uses  \code{gamlssMLpred()} and the function \code{chooseDist()} used \code{update.gamlss()}.
}
}
\value{
For the functions \code{fitDist()} and \code{fitDistPred()} a \code{gamlssML} object is return (the one which minimised the GAIC or VDEV respectively) with two extra components: 
\item{fits }{an ordered list according to the GAIC of the fitted distribution}
\item{failed}{the distributions where the \code{gamlssML)()} (or \code{gamlss()}) fits have failed}

For the function \code{chooseDist()} a matrix is returned,   with rows the different distributions and columns the different GAIC's  set by \code{k}.   

}
\references{
Rigby, R. A. and  Stasinopoulos D. M. (2005). Generalized additive models for location, scale and shape,(with discussion), 
\emph{Appl. Statist.}, \bold{54}, part 3, pp 507-554.

Rigby, R. A., Stasinopoulos, D. M.,  Heller, G. Z.,  and De Bastiani, F. (2019)
	\emph{Distributions for modeling location, scale, and shape: Using GAMLSS in R}, Chapman and Hall/CRC. An older version can be found in \url{https://www.gamlss.com/}.

Stasinopoulos D. M. Rigby R.A. (2007) Generalized additive models for location scale and shape (GAMLSS) in R.
\emph{Journal of Statistical Software}, Vol. \bold{23}, Issue 7, Dec 2007, \url{https://www.jstatsoft.org/v23/i07/}.

Stasinopoulos D. M., Rigby R.A., Heller G., Voudouris V., and De Bastiani F., (2017)
\emph{Flexible Regression and Smoothing: Using GAMLSS in R},  Chapman and Hall/CRC.  

(see also \url{https://www.gamlss.com/}).
}
\author{
Mikis Stasinopoulos, Bob Rigby, Vlasis Voudouris  and Majid Djennad.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{gamlss}}, \code{\link{gamlssML}}
}
\examples{
y <- rt(100, df=1)
m1<-fitDist(y, type="realline")
m1$fits
m1$failed
# an example of using  extra
\dontrun{
#---------------------------------------  
# Example of using the argument extra  
library(gamlss.tr)
data(tensile)
gen.trun(par=1,family="GA", type="right")
gen.trun(par=1,"LOGNO", type="right")
gen.trun(par=c(0,1),"TF", type="both")
ma<-fitDist(str, type="real0to1", trace=T,
       extra=c("GAtr", "LOGNOtr", "TFtr"), 
     data=tensile) 
ma$fits
ma$failed
#-------------------------------------
# selecting model using the prediction global deviance
# Using fitDistPred
# creating training data
y <- rt(1000, df=2)
m1 <- fitDist(y, type="realline")
m1$fits
m1$fails
# create validation data
yn <- rt(1000, df=2)
# choose distribution which fits the new data best
p1 <- fitDistPred(y, type="realline", newdata=yn)
p1$fits
p1$failed
#---------------------------------------
# using the function chooseDist()
# fitting normal distribution model
m1 <- gamlss(y~pb(x), sigma.fo=~pb(x), family=NO, data=abdom)
# choose a distribution on the real line 
# and save GAIC(k=c(2,4,6.4),  i.e. AIC, Chi-square and BIC.   
t1 <- chooseDist(m1, type="realline", parallel="snow", ncpus=4)
# the GAIC's
t1
# the distributions which failed are with NA's 
# ordering according to  BIC
getOrder(t1,3)
fm<-update(m1, family=names(getOrder(t1,3)[1]))
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{distribution}
\keyword{regression}

