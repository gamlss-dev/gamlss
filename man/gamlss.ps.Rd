\name{gamlss.ps}
\alias{gamlss.ps}
\alias{gamlss.pb}
\alias{gamlss.cy}
\alias{gamlss.pvc}
\alias{gamlss.ri}
\alias{gamlss.ridge}


\title{Support for Functions ps(), pb(), cy() and ri() }
\description{
 Those functions  are support for the functions \code{ps()} and \code{pb()}.
The functions are not intended to be called directly by users. 
}
\usage{
gamlss.ps(x, y, w, xeval = NULL, ...)
gamlss.pb(x, y, w, xeval = NULL, ...)
gamlss.ridge(x, y, w, xeval = NULL, ...)
gamlss.ri(x, y, w, xeval = NULL, ...)
gamlss.cy(x, y, w, xeval = NULL, ...)
gamlss.pvc(x, y, w, xeval = NULL, ...)
}

\arguments{
  \item{x}{ the \code{x} for function \code{gamlss.fp} is referred to the design matric of the specific parameter model (not to be used by the user)}
  \item{y}{the \code{y} for function \code{gamlss.fp} is referred to the working variable of the specific parameter model (not to be used by the user)}
  \item{w}{ the \code{w} for function \code{gamlss.fp} is referred to the iterative weight variable of the specific parameter model (not to be used by the user) }
  \item{xeval}{used in prediction}
   \item{...}{further arguments passed to or from other methods.} 
}

\value{
  
  \item{comp1 }{Description of 'comp1'}
  \item{comp2 }{Description of 'comp2'}
  ...
}
\references{ Rigby, R. A. and  Stasinopoulos D. M. (2005). Generalized additive models for location, scale and shape,(with discussion), 
\emph{Appl. Statist.}, \bold{54}, part 3, pp 507-554.

Stasinopoulos D. M., Rigby R.A. and Akantziliotou C. (2003) Instructions on how to use the GAMLSS package in R.
Accompanying documentation in the current GAMLSS  help files, (see also  \url{http://www.gamlss.org/}).  

Stasinopoulos D. M. Rigby R.A. (2007) Generalized additive models for location scale and shape (GAMLSS) in R.
\emph{Journal of Statistical Software}, Vol. \bold{23}, Issue 7, Dec 2007, \url{http://www.jstatsoft.org/v23/i07}.
}

\author{ Mikis Stasinopoulos \email{d.stasinopoulos@londonmet.ac.uk}, Bob Rigby \email{r.rigby@londonmet.ac.uk} }

\seealso{ \code{\link{gamlss}}, \code{\link{ps}},  \code{\link{pb}} }

\keyword{regression}% 
