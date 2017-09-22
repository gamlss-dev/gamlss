\name{fitDist}
\alias{fitDist}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Fits Different Parametric \code{gamlss.family} distributions to data
}
\description{
This function is using the function \code{gamlssML()} to fit all relevant parametric \code{gamlss.family} distributions to a data vector.  The final model is the one which is selected by the generalised Akaike information criterion with penalty \code{k}.

}
\usage{
fitDist(y, k = 2, 
      type = c("realAll", "realline", "realplus", "real0to1", "counts", "binom"), 
      try.gamlss = FALSE, extra = NULL, data = NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{the data vector}
  \item{k}{the penalty for the GAIC with default values \code{k=2} the standard AIC}
  \item{type}{the type of distribution to be tried see details}
  \item{try.gamlss}{if \code{gamlssML()} failed whether should try \code{gamlss} instead. This will slow up things for big data.}
  \item{extra}{whether extra distribution should be tried which are not in the \code{type} list}
  \item{data}{the data frame where \code{y} ca be found}
  \item{\dots}{for extra arguments to be passed to gamlssML() to gamlss()}
}
\details{
The following are the different \code{type} argument:
\itemize{
  \item{\code{realAll}}{ all the \code{gamlss.family} continuous distributions defined on the real line, i.e. \code{realline} plus \code{realplus}}
  \item{\code{realline}}{the \code{gamlss.family} continuous distributions : "GU", "RG" ,"LO", "NET", "TF", "PE",  "SN1",   "SN2", "SHASH", "EGB2", "JSU", "SEP1", "SEP2", "SEP3", "SEP4","ST1", "ST2", "ST3", "ST4", "ST5", "GT"}
  \item{\code{realplus}}{ the \code{gamlss.family} continuous distributions in the positive leal line: "EXP","GA","IG","LNO", "WEI3", "BCCGo", "exGAUS", "GG", "GIG", "BCTo", "BCPEo"}
  \item{\code{real0to1}}{the \code{gamlss.family} continuous distributions from 0 to 1: "BE", "BEINF", "BEINF0", "BEINF1", "BEOI", "BEZI", "GB1"}
 \item{\code{counts}}{the \code{gamlss.family} distributions for counts: "PO", "LG", "NBI", "NBII", "PIG", "DEL", "SI", "ZIP", "ZAP", "ZALG", "ZANBI", "ZIP2", "ZIPIG"}
\item{\code{binom}}{the \code{gamlss.family} distributions for binomial type data :"BI", "BB", "ZIBI", "ZIBB", "ZABI", "ZABB"}
}
}
\value{
A \code{gamlssML} object with two extra components: 
\item{fits }{an ordered list according to the GAIC of the fitted distribution}
\item{failed}{the distributions where the \code{gamlssML)()} (or \code{gamlss()}) fits have failed}
}
\references{
Rigby, R. A. and  Stasinopoulos D. M. (2005). Generalized additive models for location, scale and shape,(with discussion), 
\emph{Appl. Statist.}, \bold{54}, part 3, pp 507-554.

Stasinopoulos D. M. Rigby R.A. (2007) Generalized additive models for location scale and shape (GAMLSS) in R.
\emph{Journal of Statistical Software}, Vol. \bold{23}, Issue 7, Dec 2007, \url{http://www.jstatsoft.org/v23/i07}.

Stasinopoulos D. M., Rigby R.A., Heller G., Voudouris V., and De Bastiani F., (2017) \emph{Flexible Regression and Smoothing: Using GAMLSS in R}, Chapman and Hall/CRC. 

(see also \url{http://www.gamlss.org/}).
}
\author{
Mikis Stasinopoulos \email{mikis.stasinopoulos@gamlss.org}, Bob Rigby, Vlasis Voudouris \email{v.voiudouris@londonmet.ac.uk} and Majid Djennad \email{m.djennad.londonmet.ac.uk}
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{gamlss}},\code{\link{gamlssML}}
}
\examples{
y <- rt(100, df=1)
m1<-fitDist(y, type="realline")
m1$fits
m1$failed
# an example of using  extra
\dontrun{ 
library(gamlss.tr)
data(tensile)
gen.trun(par=1,family="GA", type="right")
gen.trun(par=1,"LOGNO", type="right")
gen.trun(par=c(0,1),"TF", type="both")
ma<-fitDist(str, type="real0to1", extra=c("GAtr", "LOGNOtr", "TFtr"), data=tensile)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{distribution}
\keyword{regression}

