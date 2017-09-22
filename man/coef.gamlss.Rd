\name{coef.gamlss}
\alias{coef.gamlss}
\title{Extract Model Coefficients in a GAMLSS fitted model}
\description{\code{coef.gamlss} is the GAMLSS specific method for the generic function \code{coef} which extracts model coefficients
     from objects returned by modelling functions.  `coefficients' is an
     alias for \code{coef}.
}
\usage{
\method{coef}{gamlss}(object, what = c("mu", "sigma", "nu", "tau"), 
                      parameter = NULL, ... )
}
\arguments{
  \item{object}{a GAMLSS fitted model}
  \item{what}{which parameter coefficient is required, default \code{what="mu"} }
   \item{parameter}{equivalent to \code{what} (more obvious name)}
  \item{\dots}{for extra arguments}
}

\value{
 Coefficients extracted from the GAMLSS model object.
}
\references{Rigby, R. A. and  Stasinopoulos D. M. (2005). Generalized additive models for location, scale and shape,(with discussion), 
\emph{Appl. Statist.}, \bold{54}, part 3, pp 507-554.

Stasinopoulos D. M. Rigby R.A. (2007) Generalized additive models for location scale and shape (GAMLSS) in R.
\emph{Journal of Statistical Software}, Vol. \bold{23}, Issue 7, Dec 2007, \url{http://www.jstatsoft.org/v23/i07}.

Stasinopoulos D. M., Rigby R.A., Heller G., Voudouris V., and De Bastiani F., (2017) \emph{Flexible Regression and Smoothing: Using GAMLSS in R}, Chapman and Hall/CRC. 

(see also \url{http://www.gamlss.org/}).
}
\author{ Mikis Stasinopoulos \email{mikis.stasinopoulos@gamlss.org}}


 
\seealso{  \code{\link{gamlss}}, \code{\link{deviance.gamlss}},  \code{\link{fitted.gamlss}} }
\examples{
data(aids)
h<-gamlss(y~poly(x,3)+qrt, family=PO, data=aids) # 
coef(h)
rm(h)
}
\keyword{regression}% 
