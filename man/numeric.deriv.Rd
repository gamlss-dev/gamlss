\name{numeric.deriv}
\alias{numeric.deriv}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{An internal GAMLSS function for numerical derivatives}
\description{
A function to calculate numerical derivatives.
}
\usage{
numeric.deriv(expr, theta, delta = NULL, 
              rho = sys.frame(sys.parent()))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{expr}{The expression to be differentiated}
  \item{theta}{A character vector}
  \item{delta}{constant for the accuracy}
  \item{rho}{environment}
}
\details{
This function is use by several GAMLSS functions  but it is not for general use since there are more reliable  function to do that in \code{R}.
}
\value{
A vector of numerical derivatives
}
\references{
Rigby, R. A. and  Stasinopoulos D. M. (2005). Generalized additive models for location, scale and shape,(with discussion),  \emph{Appl. Statist.}, \bold{54}, part 3, pp 507-554.

Rigby, R. A., Stasinopoulos, D. M.,  Heller, G. Z.,  and De Bastiani, F. (2019)
	\emph{Distributions for modeling location, scale, and shape: Using GAMLSS in R}, Chapman and Hall/CRC. An older version can be found in \url{https://www.gamlss.com/}.

Stasinopoulos D. M. Rigby R.A. (2007) Generalized additive models for location scale and shape (GAMLSS) in R.
\emph{Journal of Statistical Software}, Vol. \bold{23}, Issue 7, Dec 2007, \url{https://www.jstatsoft.org/v23/i07/}.

Stasinopoulos D. M., Rigby R.A., Heller G., Voudouris V., and De Bastiani F., (2017)
\emph{Flexible Regression and Smoothing: Using GAMLSS in R},  Chapman and Hall/CRC.  

(see also \url{https://www.gamlss.com/}).
}
\author{Mikis Stasinopoulos}
\note{
Do not use this function unless you know what you are doing
}



% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{regression}% use one of  RShowDoc("KEYWORDS")

