\name{gamlss.ps}
\alias{gamlss.ps}
\alias{gamlss.pb}
\alias{gamlss.pbo}
\alias{gamlss.cy}
\alias{gamlss.pvc}
\alias{gamlss.ri}
\alias{gamlss.pbm}
\alias{gamlss.pbz}
\alias{gamlss.pbc}
\alias{gamlss.pbp}

\title{Support for Functions  for smoothers}
\description{
 Those functions  are support for the functions  \code{pb()}, \code{pbo()},  \code{ps()},  \code{ridge()},  \code{ri()},  \code{cy()},  \code{pvc()},  and \code{pbm()}.
The functions are not intended to be called directly by users. 
}
\usage{
gamlss.pb(x, y, w, xeval = NULL, ...)
gamlss.pbo(x, y, w, xeval = NULL, ...)
gamlss.ps(x, y, w, xeval = NULL, ...)
gamlss.ri(x, y, w, xeval = NULL, ...)
gamlss.cy(x, y, w, xeval = NULL, ...)
gamlss.pvc(x, y, w, xeval = NULL, ...)
gamlss.pbm(x, y, w, xeval = NULL, ...)
gamlss.pbz(x, y, w, xeval = NULL, ...)
gamlss.pbc(x, y, w, xeval = NULL, ...)
gamlss.pbp(x, y, w, xeval = NULL, ...)
}

\arguments{
  \item{x}{ the \code{x} for function \code{gamlss.fp} is referred to the design matric of the specific parameter model (not to be used by the user)}
  \item{y}{the \code{y} for function \code{gamlss.fp} is referred to the working variable of the specific parameter model (not to be used by the user)}
  \item{w}{ the \code{w} for function \code{gamlss.fp} is referred to the iterative weight variable of the specific parameter model (not to be used by the user) }
  \item{xeval}{used in prediction}
   \item{...}{further arguments passed to or from other methods.} 
}

\value{
All function return fitted smoothers.  
}
\references{ Rigby, R. A. and  Stasinopoulos D. M. (2005). Generalized additive models for location, scale and shape,(with discussion), 
\emph{Appl. Statist.}, \bold{54}, part 3, pp 507-554.


Stasinopoulos D. M. Rigby R.A. (2007) Generalized additive models for location scale and shape (GAMLSS) in R.
\emph{Journal of Statistical Software}, Vol. \bold{23}, Issue 7, Dec 2007, \url{http://www.jstatsoft.org/v23/i07}.

Stasinopoulos D. M., Rigby R.A., Heller G., Voudouris V., and De Bastiani F., (2017) \emph{Flexible Regression and Smoothing: Using GAMLSS in R}, Chapman and Hall/CRC. 

(see also \url{http://www.gamlss.org/}).
}

\author{ Mikis Stasinopoulos \email{mikis.stasinopoulos@gamlss.org}, Bob Rigby}

\seealso{ \code{\link{gamlss}}, \code{\link{pb}},  \code{\link{ps}}, \code{\link{ri}},\code{\link{ridge}},\code{\link{cy}},\code{\link{pvc}},\code{\link{pbm}} }

\keyword{regression}% 
