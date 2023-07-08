\name{gamlssML}
\alias{gamlssML}
\alias{gamlssMLpred}


\title{Maximum Likelihood estimation of a simple GAMLSS model }
\description{
 The function \code{gamlssML()} fits  a \code{gamlss.family} distribution to single  data set  using a non linear maximisation algorithm in \code{R}.
This is relevant only when explanatory variables do not exist.

The function \code{gamlssMLpred()} is similar to \code{gamlssML()} but it saves the \emph{predictive} global deviance for the \code{newdata}.  The new data in \code{gamlssMLpred()} can be given with the arguments \code{newdata} or defining the factor  \code{rand}.  \code{rand} should be a binary factor \code{rand} splitting the original  data set into a training set (value 1)  and  a validation/test set (values 2), see
also \code{\link{gamlssVGD}}} 
\usage{
gamlssML(formula, family = NO, weights = NULL, mu.start = NULL, 
 sigma.start = NULL, nu.start = NULL, tau.start = NULL, 
 mu.fix = FALSE, sigma.fix = FALSE, nu.fix = FALSE, 
 tau.fix = FALSE, data, start.from = NULL, ...)

gamlssMLpred(response = NULL, data = NULL, family = NO, 
 rand = NULL, newdata = NULL, ...) 
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula, response}{
a vector of data requiring the fit of a \code{gamlss.family} distribution or (only for the function \code{gamlssML}) a formula,  for example, \code{y~1}, with no explanatory variables because they are ignored).  
}
  \item{family}{ \code{\link[gamlss.dist]{gamlss.family}} object, which is used to define the distribution and the link functions of the various parameters. 
   The distribution families supported by \code{gamlssML()} can be found in \code{\link[gamlss.dist]{gamlss.family}}}
  \item{weights}{ a vector of weights. 
  Here weights can be used to weight out observations (like in \code{subset}) or for a weighted likelihood analysis where  the contribution of the observations to the likelihood differs according to \code{weights}. The length of \code{weights} must be the same as the number of observations in the data. By default, the weight is set to one. To set  weights to vector say \code{w} use \code{weights=w}
}
  \item{mu.start}{ a scalar of initial values for the location parameter \code{mu} e.g.  \code{mu.start=4}}
  \item{sigma.start}{a scalar of initial values for the scale parameter \code{sigma} e.g. \code{sigma.start=1}}
  \item{nu.start}{scalar of initial values for the parameter \code{nu} e.g. \code{nu.start=3} }
  \item{tau.start}{scalar of initial values for the parameter \code{tau} e.g. \code{tau.start=3} }
  \item{mu.fix}{whether the mu parameter should be kept fixed in the fitting processes e.g. \code{mu.fix=FALSE} }
  \item{sigma.fix}{whether the sigma parameter should be kept fixed in the fitting processes e.g. \code{sigma.fix=FALSE}}
  \item{nu.fix}{whether the nu parameter should be kept fixed in the fitting processes e.g. \code{nu.fix=FALSE} }
  \item{tau.fix}{whether the tau parameter should be kept fixed in the fitting processes e.g. \code{tau.fix=FALSE}}
  \item{data}{a data frame containing the variable \code{y}, e.g. \code{data=aids}. If this is missing, the variable should be on the search list.}
 \item{start.from}{a gamlss object to start from the fitting or vector of length as many parameters in the distribution}
 \item{rand}{For \code{gamlssMLpred()} a factor with values 1 (for fitting) and 2 (for predicting).}
  \item{newdata}{The prediction data  set (validation or test).}
  \item{\dots}{for extra arguments}
}
\details{ The  function \code{gamlssML()} fits a \code{gamlss.family} distribution to a single data set is using a non linear maximisation.
 in fact it uses the internal function \code{MLE()} which is a copy of the \code{mle()} function of package \code{stat4}.
The function \code{gamlssML()} could be for large data faster than the equivalent \code{gamlss()} function which is designed for regression type of models. 

The  function \code{gamlssMLpred()} uses the  function  \code{gamlssML()}   to fit the model but then uses \code{predict.gamlssML()} to predict for new data and saves the the prediction i) deviance increments, ii) global deviance iii) residuals.   
}
\value{Returns a \code{gamlssML} object which behaves like a \code{gamlss} fitted objected}
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
Mikis Stasinopoulos \email{d.stasinopoulos@londonmet.ac.uk}, Bob Rigby, Vlasis Voudouris and Majid Djennad 
}

\seealso{
 \code{\link[gamlss.dist]{gamlss.family}}, \code{\link{gamlss}} }
\examples{
#-------- negative binomial 1000 observations
y<- rNBI(1000)
  system.time(m1<-gamlss(y~1, family=NBI))
  system.time(m1a<-gamlss(y~1, family=NBI, trace=FALSE))
system.time(m11<-gamlssML(y, family=NBI))
AIC(m1,m1a,m11, k=0)
# neg. binomial   n=10000
 y<- rNBI(10000)
 system.time(m1<-gamlss(y~1, family=NBI))
 system.time(m1a<-gamlss(y~1, family=NBI, trace=FALSE))
system.time(m11<-gamlssML(y, family=NBI))
AIC(m1,m1a,m11, k=0)
# binomial type data 
data(aep)
m1 <- gamlssML(aep$y, family=BB) # ok
m2 <- gamlssML(y, data=aep, family=BB) # ok
m3 <- gamlssML(y~1, data=aep, family=BB) # ok 
m4 <- gamlssML(aep$y~1, family=BB) # ok
AIC(m1,m2,m3,m4)
\dontrun{
#-----------------------------------------------------------
# neg. binomial   n=10000
y<- rNBI(10000)
rand <- sample(2, length(y), replace=TRUE, prob=c(0.6,0.4))
table(rand)
   Y <- subset(y, rand==1)
YVal <- subset(y, rand==2) 
length(Y)
length(YVal) 
da1 <- data.frame(y=y)
dim(da1)
da2 <- data.frame(y=Y)
dim(da2)
danew <- data.frame(y=YVal)
# using gamlssVGD to fit the models
g1 <- gamlssVGD(y~1, rand=rand, family=NBI, data=da1)
g2 <- gamlssVGD(y~1, family=NBI, data=da2, newdata=dan)
AIC(g1,g2)
VGD(g1,g2)
# using gamlssMLpred to fit the models
p1 <- gamlssMLpred(y, rand=rand, family=NBI)
p2 <- gamlssMLpred(Y, family=NBI, newdata=YVal)
# AIC and VGD should produce identical results
AIC(p1,p2,g1,g2)
VGD(p1,p2, g1,g2)
# the fitted residuals
wp(p1, ylim.all=1)
# the prediction residuals 
wp(resid=p1$residVal, ylim.all=.5)
#-------------------------------------------------------------
# chossing between distributions
p2<-gamlssMLpred(y, rand=rand, family=PO)
p3<-gamlssMLpred(y, rand=rand, family=PIG)
p4<-gamlssMLpred(y, rand=rand, family=BNB)
AIC(p1, p2, p3, p4)
VGD(p1, p2, p3, p4)
#--------------------------------------------------
}
 }

\keyword{regression}
