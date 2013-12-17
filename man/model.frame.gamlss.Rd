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
\method{model.frame}{gamlss}(formula, what = c("mu", "sigma", "nu", "tau"), ...)
\method{terms}{gamlss}(x, what = c("mu", "sigma", "nu", "tau"), ...)
\method{model.matrix}{gamlss}(object, what = c("mu", "sigma", "nu", "tau"), ...)
}
\arguments{
  \item{formula}{a gamlss object}
  \item{x}{a gamlss object}
  \item{object}{a gamlss object}
  \item{what}{for which parameter to extract the model.frame, terms or model.frame }
  \item{\dots}{for extra arguments}
}

\value{
  a model.frame, a model.matrix or terms 
}
\references{

Rigby, R. A. and  Stasinopoulos D. M. (2005). Generalized additive models for location, scale and shape,(with discussion), 
\emph{Appl. Statist.}, \bold{54}, part 3, pp 507-554.

Stasinopoulos D. M., Rigby R.A. and Akantziliotou C. (2006) Instructions on how to use the GAMLSS package in R.
Accompanying documentation in the current GAMLSS  help files, (see also  \url{http://www.gamlss.org/}). 

Stasinopoulos D. M. Rigby R.A. (2007) Generalized additive models for location scale and shape (GAMLSS) in R.
\emph{Journal of Statistical Software}, Vol. \bold{23}, Issue 7, Dec 2007, \url{http://www.jstatsoft.org/v23/i07}.
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
