\name{gamlss.random}
\alias{gamlss.random}

\title{Support for Function random()}
\description{
 This is support for the function random().
It is not intended to be called directly by users.
The function \code{gamlss.radom} is similar to the GAMLSS function \code{gamlss.ra}.
}
\usage{
gamlss.random(x, y, w)
}

\arguments{
  \item{x}{the explanatory design matrix }
  \item{y}{the response variable }
  \item{w}{iterative weights }
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

Stasinopoulos D. M., Rigby R.A. and Akantziliotou C. (2006) Instructions on how to use the GAMLSS package in R.
Accompanying documentation in the current GAMLSS  help files, (see also  \url{http://www.gamlss.org/}).

Stasinopoulos D. M. Rigby R.A. (2007) Generalized additive models for location scale and shape (GAMLSS) in R.
\emph{Journal of Statistical Software}, Vol. \bold{23}, Issue 7, Dec 2007, \url{http://www.jstatsoft.org/v23/i07}.
}
\author{Mikis Stasinopoulos, based on Trevor Hastie function \code{gam.random}}


\seealso{ \code{\link{gamlss}}, \code{\link{random}}, \code{\link{ra}}  }


\keyword{regression}% 
