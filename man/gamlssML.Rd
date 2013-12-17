\name{gamlssML}
\alias{gamlssML}

\title{Maximum Likelihood estimation of a simple GAMLSS model }
\description{
 This is a function for fitting  a \code{gamlss.family} distribution to single  data set  using a non linear maximisation algorithm in \code{R}.
This is relevant only when there are not explanatory variables. 
}
\usage{
gamlssML(y, family = NO, weights = NULL, mu.start = NULL, 
 sigma.start = NULL, nu.start = NULL, tau.start = NULL, 
 mu.fix = FALSE, sigma.fix = FALSE, nu.fix = FALSE, 
 tau.fix = FALSE, data = NULL, start.from = NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{
a vector of data requiring the fit of a \code{gamlss.family} distribution or a formula,  for example, \code{y~1} (explanatory variables are ignored).  
}
  \item{family}{ \code{\link[gamlss.dist]{gamlss.family}} object, which is used to define the distribution and the link functions of the various parameters. 
   The distribution families supported by \code{gamlssML()} can be found in \code{\link[gamlss.dist]{gamlss.family}}}
  \item{weights}{ a vector of weights. 
  Here weights can be used to weight out observations (like in \code{subset}) or for a weighted likelihood analysis where  the contribution of the observations to the likelihood differs according to \code{weights}. The length of \code{weights} must be the same as the number of observations in the data. By default, the weight is set to one. To set  weights to vector w use \code{weights=w}
}
  \item{mu.start}{ a scalar of initial values for the location parameter \code{mu} e.g.  \code{mu.start=4}}
  \item{sigma.start}{a scalar of initial values for the scale parameter \code{sigma} e.g. \code{sigma.start=1}}
  \item{nu.start}{scalar of initial values for the parameter \code{nu} e.g. \code{nu.start=3} }
  \item{tau.start}{scalar of initial values for the parameter \code{tau} e.g. \code{tau.start=3} }
  \item{mu.fix}{whether the mu parameter should be kept fixed in the fitting processes e.g. \code{mu.fix=FALSE} }
  \item{sigma.fix}{whether the sigma parameter should be kept fixed in the fitting processes e.g. \code{sigma.fix=FALSE}}
  \item{nu.fix}{whether the nu parameter should be kept fixed in the fitting processes e.g. \code{nu.fix=FALSE} }
  \item{tau.fix}{whether the tau parameter should be kept fixed in the fitting processes e.g. \code{tau.fix=FALSE}}
  \item{data}{a data frame containing the variable \code{y}. If this is missing, the variable should be on the search list.
     e.g. \code{data=aids}}
 \item{start.from}{a gamlss object to start from the fitting or vector of length as many parameters in the distribution}
  \item{\dots}{for extra arguments}
}
\details{ This function which fits a \code{gamlss.family} distribution to a single data set is using a non linear maximisation.
 in fact it uses the internal function \code{MLE()} which is a copy of the \code{mle()} function of package \code{stat4}.
The function \code{gamlssML()} could be for large data faster than the equivalent \code{gamlss()} function which is designed for regression type of models. 
}
\value{Returns a \code{gamlssML} object which behaves like a \code{gamlss} fitted objected}
\references{

Rigby, R. A. and  Stasinopoulos D. M. (2005). Generalized additive models for location, scale and shape,(with discussion), 
\emph{Appl. Statist.}, \bold{54}, part 3, pp 507-554.

Stasinopoulos D. M., Rigby R.A. and Akantziliotou C. (2006) Instructions on how to use the GAMLSS package in R.
Accompanying documentation in the current GAMLSS  help files, (see also  \url{http://www.gamlss.org/}).

Stasinopoulos D. M. Rigby R.A. (2007) Generalized additive models for location scale and shape (GAMLSS) in R.
\emph{Journal of Statistical Software}, Vol. \bold{23}, Issue 7, Dec 2007, \url{http://www.jstatsoft.org/v23/i07}.

}
\author{
Mikis Stasinopoulos \email{d.stasinopoulos@londonmet.ac.uk}, Bob Rigby \email{r.rigby@londonmet.ac.uk}, Vlasis Voudouris \email{v.voiudouris@londonmet.ac.uk} and Majid Djennad \email{m.djennad.londonmet.ac.uk}
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
m3 <- gamlssML(y~1, data=aep, family=BB) # ok but not 
#gamlssML(aep$y~1, family=BB)  
 }

\keyword{regression}
