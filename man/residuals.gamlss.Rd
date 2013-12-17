\name{residuals.gamlss}
\alias{residuals.gamlss}

\title{Extract Residuals from GAMLSS model}
\description{
 \code{residuals.gamlss} is the GAMLSS specific method for the generic function \code{residuals} which extracts the  
      residuals for a fitted model.  The abbreviated form \code{resid} is an alias for \code{residuals}.
}
\usage{
\method{residuals}{gamlss}(object, what = c("z-scores", "mu", "sigma", "nu", "tau"), 
                         type = c("simple", "weighted", "partial"), 
                         terms=NULL, ...)
}

\arguments{
  \item{object}{a GAMLSS fitted model}
  \item{what}{specify whether the standardized residuals are required, called here the "z-scores", or residuals for a specific parameter}
  \item{type}{the type of residual if residuals for a parameter are required}
  \item{terms}{if type is "partial" this specifies which term is required}
  \item{\dots}{for extra arguments}
}
\details{The "z-scores" residuals saved in a GAMLSS object are the normalized (randomized) quantile residuals (see Dunn and Smyth, 1996). 
Randomization is only needed for the discrete family distributions, see also \code{\link{rqres.plot}}. Residuals for a specific parameter can be
"simple" = (working variable - linear predictor), "weighted"= sqrt(working weights)*(working variable - linear predictor) or
 "partial"=  (working variable - linear predictor)+contribution of specific terms. 
} 
\value{
a vector or a matrix of the appropriate residuals of a GAMLSS model. Note that when weights are used in the fitting the length of the residuals can be 
different from \code{N} the length of the fitted values. Observations with weights equal to zero are not appearing in the residuals. 
Also observations with frequencies as weights will appear more than once according to their frequencies.      
}
\references{ 
Rigby, R. A. and  Stasinopoulos D. M. (2005). Generalized additive models for location, scale and shape,(with discussion), 
\emph{Appl. Statist.}, \bold{54}, part 3, pp 507-554.

Stasinopoulos D. M., Rigby R.A. and Akantziliotou C. (2006) Instructions on how to use the GAMLSS package in R.
Accompanying documentation in the current GAMLSS  help files, (see also  \url{http://www.gamlss.org/}).

Stasinopoulos D. M. Rigby R.A. (2007) Generalized additive models for location scale and shape (GAMLSS) in R.
\emph{Journal of Statistical Software}, Vol. \bold{23}, Issue 7, Dec 2007, \url{http://www.jstatsoft.org/v23/i07}.
}


\author{Mikis Stasinopoulos \email{d.stasinopoulos@londonmet.ac.uk}, Bob Rigby \email{r.rigby@londonmet.ac.uk}}
\note{The "weighted" residuals of a specified parameter can be zero and one if the square of first derivative 
     have been used in the fitting of this parameter}
\seealso{ \code{\link{print.gamlss}}, \code{\link{summary.gamlss}},  \code{\link{fitted.gamlss}},  \code{\link{coef.gamlss}},  
          \code{\link{residuals.gamlss}}, \code{\link{update.gamlss}},  \code{\link{plot.gamlss}}, \code{\link{deviance.gamlss}}, \code{\link{formula.gamlss}} 
           }
\examples{
data(aids)
h<-gamlss(y~poly(x,3)+qrt, family=NBI, data=aids) # 
plot(aids$x,resid(h))
plot(aids$x,resid(h,"sigma") )
rm(h)
}
\keyword{regression}% 
