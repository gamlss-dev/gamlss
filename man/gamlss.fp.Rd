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
  \item{x}{the \code{x} for function \code{gamlss.fp} is referred to the design matrix of the specific parameter model (not to be used by the user)}
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

Rigby, R. A., Stasinopoulos, D. M.,  Heller, G. Z.,  and De Bastiani, F. (2019)
	\emph{Distributions for modeling location, scale, and shape: Using GAMLSS in R}, Chapman and Hall/CRC. An older version can be found in \url{https://www.gamlss.com/}.

Stasinopoulos D. M. Rigby R.A. (2007) Generalized additive models for location scale and shape (GAMLSS) in R.
\emph{Journal of Statistical Software}, Vol. \bold{23}, Issue 7, Dec 2007, \url{https://www.jstatsoft.org/v23/i07/}.

Stasinopoulos D. M., Rigby R.A., Heller G., Voudouris V., and De Bastiani F., (2017)
\emph{Flexible Regression and Smoothing: Using GAMLSS in R},  Chapman and Hall/CRC.  

(see also \url{https://www.gamlss.com/}).
}

\author{ Mikis Stasinopoulos, Bob Rigby}

\seealso{ \code{\link{gamlss}}, \code{\link{fp}} }

\keyword{regression}% 
