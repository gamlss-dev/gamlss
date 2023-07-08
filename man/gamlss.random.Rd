\name{gamlss.random}
\alias{gamlss.random}
\alias{gamlss.re}

\title{Support for Functions random() and re()}

\description{
 This is support for the functions \code{random()} and \code{re()} respectively.
It is not intended to be called directly by users.
.
}
\usage{
gamlss.random(x, y, w, xeval = NULL, ...)
gamlss.re(x, y, w, xeval = NULL, ...)
}

\arguments{
  \item{x}{the explanatory design matrix }
  \item{y}{the response variable }
  \item{w}{iterative weights }
  \item{xeval}{it used internaly for prediction}
  \item{\dots}{for extra arguments}
}

\value{Returns a list with
 
  \item{y}{the fitted values}
  \item{residuals}{the residuals}
  \item{var}{the variance of the fitted values}
  \item{lambda}{the final lambda, the smoothing parameter}
  \item{coefSmo}{with value NULL}
}
\references{ Chambers, J. M. and Hastie, T. J. (1991). \emph{Statistical Models in S}, Chapman and Hall, London. 

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
\author{Mikis Stasinopoulos, based on Trevor Hastie function \code{gam.random}}


\seealso{ \code{\link{gamlss}}, \code{\link{random}}}


\keyword{regression}% 
