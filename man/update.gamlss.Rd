\name{update.gamlss}
\alias{update.gamlss}

\title{ Update and Re-fit a GAMLSS Model}
\description{\code{update.gamlss} is the GAMLSS specific method for the generic function \code{update} 
which updates and (by default) refits a GAMLSS model.
}
\usage{
\method{update}{gamlss}(object, formula., ..., what = c("mu", "sigma", "nu", "tau", "All"), 
              evaluate = TRUE)
}

\arguments{
  \item{object}{a GAMLSS fitted model}
  \item{formula.}{the formula to update}
  \item{\dots}{for updating argument in \code{gamlss()}}
  \item{what}{the parameter in which the formula needs updating for example "mu", "sigma", "nu" "tau" or "All".
   If "All" all the formulae are updated. Note that the \code{what} argument has an effect only if only if the 
   argument \code{formula.} is set}
  \item{evaluate}{whether to evaluate the call or not}
}

\value{
Returns a GAMLSS call or fitted object.  
}
\references{ Rigby, R. A. and  Stasinopoulos D. M. (2005). Generalized additive models for location, scale and shape,(with discussion), 
\emph{Appl. Statist.}, \bold{54}, part 3, pp 507-554.

Stasinopoulos D. M., Rigby R.A. and Akantziliotou C. (2006) Instructions on how to use the GAMLSS package in R.
Accompanying documentation in the current GAMLSS  help files, (see also  \url{http://www.gamlss.org/}).

Stasinopoulos D. M. Rigby R.A. (2007) Generalized additive models for location scale and shape (GAMLSS) in R.
\emph{Journal of Statistical Software}, Vol. \bold{23}, Issue 7, Dec 2007, \url{http://www.jstatsoft.org/v23/i07}.
 
}
\author{ Mikis Stasinopoulos \email{d.stasinopoulos@londonmet.ac.uk}, Bob Rigby \email{r.rigby@londonmet.ac.uk}  }
\seealso{ \code{\link{print.gamlss}}, \code{\link{summary.gamlss}},  \code{\link{fitted.gamlss}},  \code{\link{coef.gamlss}},  
          \code{\link{residuals.gamlss}},  \code{\link{plot.gamlss}}, \code{\link{deviance.gamlss}}, \code{\link{formula.gamlss}}  }
\examples{
data(aids)
# fit a poisson model
h.po <-gamlss(y~pb(x)+qrt, family=PO, data=aids) 
# update with a negative binomial
h.nb <-update(h.po, family=NBI) 
# update the smoothing 
h.nb1 <-update(h.nb,~cs(x,8)+qrt) 
# remove qrt
h.nb2 <-update(h.nb1,~.-qrt)
# put back qrt take log of y and fit a normal distribution 
h.nb3 <-update(h.nb1,log(.)~.+qrt, family=NO) 
# verify that it is the same 
h.no<-gamlss(log(y)~cs(x,8)+qrt,data=aids )
}
\keyword{regression}% 
