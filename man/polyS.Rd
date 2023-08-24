\name{polyS}
\alias{polyS}
\alias{poly.matrix}

\title{Auxiliary support for the GAMLSS}
\description{
These two functions are similar to the \code{poly} and \code{polym} in R. 
Are needed for the \code{gamlss.lo} function of GAMLSS and should not be used on their own.  
}
\usage{
polyS(x, ...)
poly.matrix(m, degree = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{a variable}
  \item{m}{a variable }
  \item{degree}{the degree of the polynomial}
  \item{\dots}{for extra arguments}
  
}

\value{
  Returns a matrix of orthogonal polynomials
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
\author{Mikis Stasinopoulos}


\section{Warning }{Not be use by the user}  
\seealso{ \code{\link{gamlss}}, \code{\link{gamlss.lo}} }
\keyword{regression}% 
