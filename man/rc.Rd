\name{rc}
\alias{rc}
\title{Specify Random Coefficients In A GAMLSS Formula  }
\description{Fits Random Coefficients In A GAMLSS Model.
}
\usage{
rc(formula, lambda = NULL)
}

\arguments{
  \item{formula}{a model formula to specify the explanatory variable and the subject grouping factor e.g. \code{formula=~x1|Subjects} will fit 
                 random constant and slopes for x1 for each level of Subjects  }
  \item{lambda}{a matrix specifying the variance-covariance matrix for the random coefficients e.g. \code{lambda=diag{2,5}} }
}
\details{ This is an experimental function and it is likely to change in the near future.
 For a chosen fixed variance-covariance matrix for the random coefficients, this function finds the MAP estimates for the constant and slope random effects for the different subjects
 }
\value{
The variables in the left side of the formula are  returned with class "smooth", and  attribute named "call" which is to be evaluated 
in the backfitting  \code{additive.fit()}  called by \code{gamlss()}
}
\references{
Rigby, R. A. and  Stasinopoulos D. M. (2005). Generalized additive models for location, scale and shape,(with discussion), 
\emph{Appl. Statist.}, \bold{54}, part 3, pp 507-554.

Stasinopoulos D. M., Rigby R.A. and Akantziliotou C. (2006) Instructions on how to use the GAMLSS package in R.
Accompanying documentation in the current GAMLSS  help files, (see also  \url{http://www.gamlss.org/})

Stasinopoulos D. M. Rigby R.A. (2007) Generalized additive models for location scale and shape (GAMLSS) in R.
\emph{Journal of Statistical Software}, Vol. \bold{23}, Issue 7, Dec 2007, \url{http://www.jstatsoft.org/v23/i07}.
.
}
\author{Mikis Stasinopoulos \email{d.stasinopoulos@londonmet.ac.uk}, Bob Rigby \email{r.rigby@londonmet.ac.uk} }


\section{Warning }{This is experimental} 
\seealso{ \code{\link{random}}, \code{\link{ra}}, \code{\link{gamlss}} }

\keyword{regression}% 
