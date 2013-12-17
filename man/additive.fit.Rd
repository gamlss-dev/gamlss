\name{additive.fit}
\alias{additive.fit}

\title{Implementing Backfitting in GAMLSS}
\description{
This function is not to be used on its own. It is used for backfitting in the GAMLSS fitting algorithms and it is based on 
the equivalent function written by Trevor Hastie in the gam() S-plus implementation, (Chambers and Hastie, 1991). 
}
\usage{
additive.fit(x, y, w, s, who, smooth.frame, maxit = 30, tol = 0.001, 
             trace = FALSE, se = TRUE, ...)
}

\arguments{
  \item{x}{the linear part of the explanatory variables}
  \item{y}{the response variable}
  \item{w}{the weights}
  \item{s}{the matrix containing the smoothers }
  \item{who}{the current smoothers}
  \item{smooth.frame}{the data frame used for the smoothers}
  \item{maxit}{maximum number of iterations in the backfitting}
  \item{tol}{the tolerance level for the backfitting }
  \item{trace}{whether to trace the backfitting algorithm}
  \item{se}{whether standard errors are required}
  \item{\dots}{for extra arguments}
}
\details{
  This function should not be used on its own 
}
\value{
  Returns a list with the linear fit plus the smothers 
}
\references{
Chambers, J. M. and Hastie, T. J. (1991). \emph{Statistical Models in S}, Chapman and Hall, London. 

Rigby, R. A. and  Stasinopoulos D. M. (2005). Generalized additive models for location, scale and shape,(with discussion), 
\emph{Appl. Statist.}, \bold{54}, part 3, pp 507-554.

Stasinopoulos D. M., Rigby R.A. and Akantziliotou C. (2006) Instructions on how to use the GAMLSS package in R.
Accompanying documentation in the current GAMLSS  help files, (see also  \url{http://www.gamlss.org/}).

Stasinopoulos D. M. Rigby R.A. (2007) Generalized additive models for location scale and shape (GAMLSS) in R.
\emph{Journal of Statistical Software}, Vol. \bold{23}, Issue 7, Dec 2007, \url{http://www.jstatsoft.org/v23/i07}.
}
\author{Mikis Stasinopoulos}

\seealso{ \code{\link{gamlss}} }

\keyword{regression}% 
