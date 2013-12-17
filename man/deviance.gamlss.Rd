\name{deviance.gamlss}
\alias{deviance.gamlss}

\title{ Global Deviance of a GAMLSS model }
\description{
 Returns the global, -2*log(likelihood), or the penalized, -2*log(likelihood)+ penalties, deviance of a fitted GAMLSS model object.
}
\usage{
\method{deviance}{gamlss}(object,  what = c("G", "P"), ...)
}
\arguments{
  \item{object}{a GAMLSS fitted model}
  \item{what}{put "G" for Global or "P" for Penalized deviance}
  \item{\dots}{for extra arguments}
}
\details{
\code{deviance} is a generic function which can be used to extract deviances
for fitted models.  \code{deviance.gamlss} is the method for a GAMLSS object. 
}

\value{
 The value of the global or the penalized deviance extracted from a GAMLSS object. 
}

\references{Rigby, R. A. and  Stasinopoulos D. M. (2005). Generalized additive models for location, scale and shape,(with discussion), 
\emph{Appl. Statist.}, \bold{54}, part 3, pp 507-554.

Stasinopoulos D. M., Rigby R.A. and Akantziliotou C. (2006) Instructions on how to use the GAMLSS package in R.
Accompanying documentation in the current GAMLSS  help files, (see also  \url{http://www.gamlss.org/}).   

Stasinopoulos D. M. Rigby R.A. (2007) Generalized additive models for location scale and shape (GAMLSS) in R.
\emph{Journal of Statistical Software}, Vol. \bold{23}, Issue 7, Dec 2007, \url{http://www.jstatsoft.org/v23/i07}.
}

\author{Mikis Stasinopoulos \email{d.stasinopoulos@londonmet.ac.uk} }

\seealso{ \code{\link[gamlss.dist]{gamlss.family}}, \code{\link{coef.gamlss}},  \code{\link{fitted.gamlss}} }
\examples{
data(aids)
h<-gamlss(y~poly(x,3)+qrt, family=PO, data=aids) # 
deviance(h)
rm(h)
}
\keyword{regression}% 
