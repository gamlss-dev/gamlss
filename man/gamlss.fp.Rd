\name{gamlss.fp}
\alias{gamlss.pp}
\alias{gamlss.fp}

\title{ Support for Function fp() }
\description{
Those are support for the functions \code{fp()} and \code{pp}.
It is not intended to be called directly by users. 
}
\usage{
gamlss.fp(x, y, w, npoly = 2, xeval = NULL)
gamlss.pp(x, y, w)
}

\arguments{
  \item{x}{the \code{x} for function \code{gamlss.fp} is referred to the design matric of the specific parameter model (not to be used by the user)}
  \item{y}{the \code{y} for function \code{gamlss.fp} is referred to the working variable of the specific parameter model (not to be used by the user)}
  \item{w}{the \code{w} for function \code{gamlss.fp} is referred to the iterative weight variable of the specific parameter model (not to be used by the user) }
  \item{npoly}{a positive indicating how many fractional polynomials should be considered in the fit. Can take the values 1, 2 or 3 with 2 as default }
  \item{xeval}{used in prediction } 
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

\seealso{ \code{\link{gamlss}}, \code{\link{fp}} }

\keyword{regression}% 
