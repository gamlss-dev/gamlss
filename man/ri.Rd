\name{ri}
\alias{ri}
\title{Specify ridge  or lasso  Regression within a  GAMLSS Formula}
\description{
The function \code{ri()} allow the user to fit a ridge regression within GAMLSS. 
It allows the coefficients of a set of explanatory variables to be shrunk towards zero.
The amount of shrinking depends either on lambda, or on the equivalent degrees of freedom (df). The type of shrinking depends on the argument \code{Lp} see example.}
\usage{
ri(X = NULL, x.vars = NULL, df = NULL, lambda = NULL, 
   method = c("ML", "GAIC"), order = 0, start = 10, Lp = 2, 
   kappa = 1e-05, iter = 100, c.crit = 1e-06, k = 2)
}
\arguments{
  \item{X}{A matrix of explanatory variables \code{X} which is standardised (mean=0, sd=1)  automatically. Note that in order to get predictions you should use the option \code{x.vars}}
    \item{x.vars}{which variables from the \code{data.frame} declared in \code{data} needs to be included. This is a way to fit the model if predictions are required.}
  \item{df}{
the effective degrees of freedom \code{df} 
}
  \item{lambda}{
the smoothing parameter \code{lambda}
}
  \item{method}{
which method is used for the estimation of the smoothing parameter, `ML' or `GAIC' are allowed. 
}
  \item{order}{
the \code{order} of the difference applied to the coefficients with default zero. (Do not change this unless there is some ordering in the explanatory variables).)
}
  \item{start}{
starting value for lambda if it estimated using `ML' or `GAIC'
}
  \item{Lp}{
The type of penalty required, \code{Lp=2} a proper ridge regression is the default. Use \code{Lp=1} for lasso and different values for different penalties.  
}
  \item{kappa}{a regulation parameters used for the weights in the penalties.}
  \item{iter}{the number of internal iteration allowed see details.}
  \item{c.crit}{
\code{c.crit} is the convergent criterion 
}
  \item{k}{
\code{k} is the penalty if `GAIC' method is used. 
}
}
\details{
This implementation of ridge and related regressions is based on an idea of Paul Eilers which used weights in the penalty matrix. The type of weights are defined by the argument \code{Lp}. \code{Lp=2} is the standard ridge regression,  \code{Lp=1} fits a lasso regression while  \code{Lp=0} allows a "best subset"" regression see Hastie et al (2009) page 71.
}
\value{
x is returned with class "smooth", with an attribute named "call" which is to be evaluated in the backfitting  \code{additive.fit()} called by \code{gamlss()}
}
\references{

Rigby, R. A. and  Stasinopoulos D. M. (2005). Generalized additive models for location, scale and shape,(with discussion), 
\emph{Appl. Statist.}, \bold{54}, part 3, pp 507-554.

Rigby, R. and Stasinopoulos, D. M (2013) Automatic smoothing parameter selection in GAMLSS with an application to centile estimation, \emph{Statistical methods in medical research}.

Rigby, R. A., Stasinopoulos, D. M.,  Heller, G. Z.,  and De Bastiani, F. (2019)
	\emph{Distributions for modeling location, scale, and shape: Using GAMLSS in R}, Chapman and Hall/CRC. An older version can be found in \url{https://www.gamlss.com/}.

Stasinopoulos D. M. Rigby R.A. (2007) Generalized additive models for location scale and shape (GAMLSS) in R.
\emph{Journal of Statistical Software}, Vol. \bold{23}, Issue 7, Dec 2007, \url{https://www.jstatsoft.org/v23/i07/}.

Stasinopoulos D. M., Rigby R.A., Heller G., Voudouris V., and De Bastiani F., (2017)
\emph{Flexible Regression and Smoothing: Using GAMLSS in R},  Chapman and Hall/CRC.  

(see also \url{https://www.gamlss.com/}).
}
\author{
 Mikis Stasinopoulos, Bob Rigby and Paul Eilers
}
\seealso{\code{\link{gamlss}}}
\examples{
# USAIR DATA
# standarise data 1-------------------------------------------------------------
# ridge
m1<- gamlss(y~ri(x.vars=c("x1","x2","x3","x4","x5","x6")), 
            data=usair)
# lasso
m2<- gamlss(y~ri(x.vars=c("x1","x2","x3","x4","x5","x6"), Lp=1), 
     data=usair)
# best subset
m3<- gamlss(y~ri(x.vars=c("x1","x2","x3","x4","x5","x6"), Lp=0), 
     data=usair)
#--------  plotting the coefficients
op <- par(mfrow=c(3,1))
plot(getSmo(m1)) #
plot(getSmo(m2))
plot(getSmo(m3))
par(op)
}
\keyword{regression}% 