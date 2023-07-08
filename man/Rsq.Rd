\name{Rsq}
\alias{Rsq}

\title{Generalised  (Pseudo) R-squared for GAMLSS models}
\description{
This function gives the generalised R-squared of Nagelkerke (1991) for a GAMLSS model.
}
\usage{
Rsq(object, type = c("Cox Snell","Cragg Uhler","both"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{a GAMLSS object}
  \item{type}{which definition of R squared.  Can be the "Cox Snell" or the Nagelkerke, "Cragg Uhler" or "both". }
}
\details{
The \code{Rsq()} function  uses the definition for R-squared: 
\deqn{R^2=1- \left(\frac{L(0)}{L(\hat{\theta})}\right)^(2/n)}{R-squared =1-(L(0)/L(fitted))^(2/n)}
 where \eqn{L(0)} is the null model (only a constant is  fitted to all parameters) and 
 \eqn{L(\hat{\theta})}{L(fitted)} is the current fitted model. This definition sometimes is referred to as the Cox & Snell R-squared.   The Nagelkerke /Cragg & Uhler's definition divides the above with
 \deqn{1- L(0)^(2/n)}{1-L(0)^(2/n)}
}
\value{The  \code{Rsq()} produces a single value if type="Cox Snell" or "Cragg Uhler" and a list if type="both".}
\references{
Nagelkerke, N. J. (1991). A note on a general definition of the coefficient of determination. Biometrika, 78(3), 691-692.

Rigby, R. A., Stasinopoulos, D. M.,  Heller, G. Z.,  and De Bastiani, F. (2019)
	\emph{Distributions for modeling location, scale, and shape: Using GAMLSS in R}, Chapman and Hall/CRC. An older version can be found in \url{https://www.gamlss.com/}.

Stasinopoulos D. M. Rigby R.A. (2007) Generalized additive models for location scale and shape (GAMLSS) in R.
\emph{Journal of Statistical Software}, Vol. \bold{23}, Issue 7, Dec 2007, \url{https://www.jstatsoft.org/v23/i07/}.

Stasinopoulos D. M., Rigby R.A., Heller G., Voudouris V., and De Bastiani F., (2017)
\emph{Flexible Regression and Smoothing: Using GAMLSS in R},  Chapman and Hall/CRC.  

(see also \url{https://www.gamlss.com/}).
}

\author{Mikis Stasinopoulos \email{d.stasinopoulos@londonmet.ac.uk} }
\note{The null model is fitted using the function gamlssML() which can create warning messages}

\seealso{\code{\link{GAIC}}}
\examples{
data(aids)
m1 <- gamlss(y~x+qrt, data=aids, family=NBI)
Rsq(m1)
Rsq(m1, type="both")
}

\keyword{regression}

