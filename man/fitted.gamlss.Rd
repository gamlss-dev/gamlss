\name{fitted.gamlss}
\alias{fitted.gamlss}
\alias{fv}

\title{Extract Fitted Values For A GAMLSS Model }
\description{\code{fitted.gamlss} is the GAMLSS specific method for the generic function 
 \code{fitted} which extracts fitted values for a specified parameter from  a GAMLSS objects.  \code{fitted.values} is an
     alias for it. The function \code{fv()} is similar to \code{fitted.gamlls()} but allows the argument \code{what} not to be character
}
\usage{
\method{fitted}{gamlss}(object, what = c("mu", "sigma", "nu", "tau"),
                parameter= NULL, ...)
fv(obj, what = c("mu", "sigma", "nu", "tau"), parameter= NULL, ... )
}

\arguments{
  \item{object}{a GAMLSS fitted model}
  \item{obj}{a GAMLSS fitted model}
  \item{what}{which parameter fitted values  are required, default \code{what="mu"} }
  \item{parameter}{equivalent to \code{what}}
  \item{\dots}{for extra arguments}
}

\value{
  Fitted values extracted from the GAMLSS object for the given parameter.
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

\author{Mikis Stasinopoulos \email{d.stasinopoulos@londonmet.ac.uk} }


\seealso{ \code{\link{print.gamlss}}, \code{\link{summary.gamlss}},  \code{\link{fitted.gamlss}},  \code{\link{coef.gamlss}},  
          \code{\link{residuals.gamlss}}, \code{\link{update.gamlss}},  \code{\link{plot.gamlss}}, \code{\link{deviance.gamlss}}, 
          \code{\link{formula.gamlss}} }
\examples{
data(aids)
h<-gamlss(y~poly(x,3)+qrt, family=PO, data=aids) # 
fitted(h)
rm(h)
}
\keyword{regression}% 
