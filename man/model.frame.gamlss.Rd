\name{model.frame.gamlss}
\alias{model.frame.gamlss}
\alias{terms.gamlss}
\alias{model.matrix.gamlss}

\title{Extract a model.frame, a model matrix or terms from a GAMLSS object for a given distributional parameter}
\description{
  \code{model.frame.gamlss}, \code{model.matrix.gamlss} and \code{terms.gamlss} are the gamlss versions of the  generic functions
  \code{model.frame}, \code{model.matrix} and \code{terms} respectively. 
}
\usage{
\method{model.frame}{gamlss}(formula, what = c("mu", "sigma", "nu", "tau"),
                     parameter= NULL, ...)
\method{terms}{gamlss}(x, what = c("mu", "sigma", "nu", "tau"), 
                     parameter= NULL,  ...)
\method{model.matrix}{gamlss}(object, what = c("mu", "sigma", "nu", "tau"), 
                      parameter= NULL, ...)
}
\arguments{
  \item{formula}{a gamlss object}
  \item{x}{a gamlss object}
  \item{object}{a gamlss object}
  \item{what}{for which parameter to extract the model.frame, terms or model.frame }
  \item{parameter}{equivalent to \code{what}}
  \item{\dots}{for extra arguments}
}

\value{
  a model.frame, a model.matrix or terms 
}
\references{

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
\author{Mikis Stasinopoulos }
\seealso{ \code{\link{gamlss}} }
\examples{
data(aids)
mod<-gamlss(y~poly(x,3)+qrt, family=PO, data=aids) # 
model.frame(mod)
model.matrix(mod)
terms(mod, "mu")
rm(mod)
}
\keyword{regression}% 
