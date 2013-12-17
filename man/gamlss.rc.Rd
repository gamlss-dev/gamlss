\name{gamlss.rc}
\alias{gamlss.rc}

\title{ Support for Function rc() }
\description{
This is support for the function rc().
It is not intended to be called directly by users. 
}
\usage{
gamlss.rc(x, y, w, df = NULL)
}

\arguments{
  \item{x}{the explanatory design matrix}
  \item{y}{the response variable}
  \item{w}{iterative weights }
  \item{df}{effective degrees of freedom }
}

\value{
  Returns a list with
  \item{fitted.values}{fitted}
  \item{residuals}{residuals}
  \item{var}{}
  \item{nl.df}{the trace of the smoothing matrix}
  \item{lambda}{the value of the smoothing parameter}
  \item{coefSmo}{the coefficients from the smoothing fit}
  \item{varcoeff}{the variance of the coefficients}
}
\references{ Rigby, R. A. and  Stasinopoulos D. M. (2005). Generalized additive models for location, scale and shape,(with discussion), 
\emph{Appl. Statist.}, \bold{54}, part 3, pp 507-554.

Stasinopoulos D. M., Rigby R.A. and Akantziliotou C. (2006) Instructions on how to use the GAMLSS package in R.
Accompanying documentation in the current GAMLSS  help files, (see also  \url{http://www.gamlss.org/}). 

Stasinopoulos D. M. Rigby R.A. (2007) Generalized additive models for location scale and shape (GAMLSS) in R.
\emph{Journal of Statistical Software}, Vol. \bold{23}, Issue 7, Dec 2007, \url{http://www.jstatsoft.org/v23/i07}.
}

\author{ Mikis Stasinopoulos \email{d.stasinopoulos@londonmet.ac.uk}, Bob Rigby \email{r.rigby@londonmet.ac.uk} }

\seealso{ \code{\link{gamlss}}, \code{\link{rc}} }

\keyword{regression}% 
