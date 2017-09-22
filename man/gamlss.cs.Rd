\name{gamlss.cs}
\alias{gamlss.cs}


\title{Support for Function cs() and scs()}
\description{ This is support for the functions cs(), and scs().
It is not intended to be called directly by users. The function \code{gamlss.cs} is using  the R function \code{smooth.spline}
}
\usage{
gamlss.cs(x, y, w, df = NULL, spar = NULL, xeval = NULL, ...)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{x}{the design matrix}
  \item{y}{the response variable}
  \item{w}{ prior weights}
  \item{df}{ effective degrees of freedom }
  \item{spar}{spar the smoothing parameter}
  \item{xeval}{used in prediction } 
   \item{\dots}{for extra arguments}
}

\value{
  Returns a class "smooth.spline" object with 
  \item{residuals}{The residuals of the fit}
  \item{fitted.values}{The smoothing values }
  \item{var}{the variance for the fitted smoother}
  \item{lambda}{the final value for spar}
  \item{nl.df}{the smoothing degrees of freedom excluding the constant and linear terms, i.e. (df-2)}
  \item{coefSmo}{this is a list containing among others the knots and the coefficients}
}
\author{Mikis Stasinopoulos \email{mikis.stasinopoulos@gamlss.org}, Bob Rigby}

\seealso{ \code{\link{gamlss}}, \code{\link{cs}} }

\keyword{regression}% 
