\name{lo}
\alias{lo}
\alias{lo.control}
\title{ Specify a loess fit in a GAMLSS formula}
\description{
Allows the user to specify a loess fit in a GAMLSS model. This function is similar to the \code{lo} function in the \code{gam} implementation of package \code{gam} see Chambers and Hastie (1991). 
}
\usage{
lo(formula, control = lo.control(...), ...)
lo.control(span = 0.75, enp.target = NULL, 
          degree = 2, parametric = FALSE, drop.square = FALSE, 
          normalize = TRUE, family = c("gaussian", "symmetric"), 
          method = c("loess", "model.frame"), 
          surface = c("interpolate", "direct"), 
          statistics = c("approximate", "exact"), 
          trace.hat = c("exact", "approximate"), 
          cell = 0.2, iterations = 4, ...)
}

\arguments{
   \item{formula}{a formula specifying the explanatory variables }
  \item{control}{a control to be passed to the \code{loess} function }
  \item{\dots}{extra arguments}
  \item{span}{the number of observations in a neighbourhood. This is the smoothing parameter for a loess fit.}
  \item{enp.target}{an alternative way to specify span, as the approximate equivalent number degrees of freedom to be used. See also the help file of the \R{} function \code{loess}. For consistency with the older version of \code{lo} the effective degrees of freedom \code{df} can be also specified instead of span, e.g. \code{df=5}}
  \item{degree}{the degree of local polynomial; can be 1 or 2. See also the help file of \code{loess}}
  \item{parametric}{should any terms be fitted globally rather than locally? See the help file of \code{loess}}
  \item{drop.square}{for fits with more than one predictor and degree=2, should the quadratic term be dropped for particular predictors?. See also help file of \code{loess}}
  \item{normalize}{should the predictors be normalized to a common scale if there is more than one? See the help file of \code{loess}}
  \item{family}{if "gaussian" fitting is by least-squares, and if "symmetric" a re-descending M estimator is used with Tukey's biweight function. See the help file of \code{loess}}
  \item{method}{fit the model or just extract the model frame. See the help file of \code{loess}}
  \item{surface}{should the fitted surface be computed exactly or via interpolation from a kd tree? See also
  the help file of \code{loess.control}}
  \item{statistics}{should the statistics be computed exactly or approximately? See the help file of \code{loess.control}}
  \item{trace.hat}{should the trace of the smoother matrix be computed exactly or approximately? See the help file of \code{loess.control}}
  \item{cell}{if interpolation is used this controls the accuracy of the approximation via the maximum number of points in a cell in the kd tree. See the help file of \code{loess.control}}
  \item{iterations}{the number of iterations used in robust fitting. See the help file of \code{loess.control}}
}
\details{
 Note that \code{lo} itself does no smoothing; it simply sets things up for the function \code{gamlss.lo()} which is used by the backfitting function  \code{gamlss.add()}.
}
\value{
 a loess object is returned.  
}
\references{

Chambers, J. M. and Hastie, T. J. (1991). \emph{Statistical Models in S}, Chapman and Hall, London. 

Rigby, R. A. and  Stasinopoulos D. M. (2005). Generalized additive models for location, scale and shape,(with discussion), 
\emph{Appl. Statist.}, \bold{54}, part 3, pp 507-554.

Stasinopoulos D. M., Rigby R.A. and Akantziliotou C. (2006) Instructions on how to use the GAMLSS package in R.
Accompanying documentation in the current GAMLSS  help files, (see also  \url{http://www.gamlss.org/}).

Stasinopoulos D. M. Rigby R.A. (2007) Generalized additive models for location scale and shape (GAMLSS) in R.
\emph{Journal of Statistical Software}, Vol. \bold{23}, Issue 7, Dec 2007, \url{http://www.jstatsoft.org/v23/i07}.
}

\author{Mikis Stasinopoulos \email{d.stasinopoulos@londonmet.ac.uk}, Bob Rigby \email{r.rigby@londonmet.ac.uk}, (The original \code{lo()} function was based on the Trevor Hastie's S-plus \code{lo()} function. See also the documentation of the \code{loess} function for the authorship of the function.}
\note{Note that \code{lo} itself does no smoothing; it simply sets things up for \code{gamlss.lo()} to do the backfitting.}

\section{Warning }{In this version the first argument is a formula NOT a list as in the previous one} 

\seealso{ \code{\link{cs}}, \code{\link{random}},  }

\examples{
data(aids)
attach(aids)
# fitting a loess curve with span=0.4 plus the a quarterly  effect 
aids1<-gamlss(y~lo(~x,span=0.4)+qrt,data=aids,family=PO) # 
plot(x,y)
lines(x,fitted(aids1))
rm(aids1)
detach(aids)
}
\keyword{regression}% 
