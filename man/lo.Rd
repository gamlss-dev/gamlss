\name{lo}
\alias{lo}
\alias{lo.control}
\alias{vis.lo}

\title{ Specify a loess fit in a GAMLSS formula}
\description{
Allows the user to specify a \code{loess} fit within a GAMLSS model. This function is similar to the \code{lo} function in the \code{gam} implementation of package \code{gam} see Chambers and Hastie (1991).

The function \code{vis.lo()} allows plotting the results.
}
\usage{
lo(formula, control = lo.control(...), ...)
lo.control(span = 0.75, enp.target = NULL, 
          degree = 2, parametric = FALSE, drop.square = FALSE, 
          normalize = TRUE, family = c("gaussian", "symmetric"), 
          method = c("loess", "model.frame"), 
          surface = c("interpolate", "direct"), 
          statistics = c("approximate", "exact", "none"), 
          trace.hat = c("exact", "approximate"), 
          cell = 0.2, iterations = 4,iterTrace = FALSE, ...)
vis.lo(obj, se=-1, rug = FALSE, partial.resid = FALSE, 
          col.term = "darkred",  col.shaded = "gray", 
          col.res = "lightblue", col.rug = "gray", lwd.term = 1.5, 
          cex.res = 1, pch.res = par("pch"), 
          type = c("persp", "contour"), col.surface = "gray", 
          nlevels = 30, n.grid = 30, image = TRUE, ...)
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
  \item{family}{if \code{"gaussian"} fitting is by least-squares, and if "symmetric" a re-descending M estimator is used with Tukey's biweight function. See the help file of \code{loess}}
  \item{method}{fit the model or just extract the model frame. See the help file of \code{loess}}
  \item{surface}{should the fitted surface be computed exactly or via interpolation from a kd tree? See also
  the help file of \code{loess.control}}
  \item{statistics}{should the statistics be computed exactly or approximately? See the help file of \code{loess.control}}
  \item{trace.hat}{should the trace of the smoother matrix be computed exactly or approximately? See the help file of \code{loess.control}}
  \item{cell}{if interpolation is used this controls the accuracy of the approximation via the maximum number of points in a cell in the kd tree. See the help file of \code{loess.control}}
  \item{iterations}{the number of iterations used in robust fitting. See the help file of \code{loess.control}}
\item{iterTrace}{logical (or integer) determining if tracing information during the robust iterations (iterations>= 2) is produced. See the help file of \code{loess.control}}
%\item{se}{logical, indicating whether the estimation of the standard errors for the smooth function should be estimated within the gamlss backfitting or not. NOTE that \code{se=TRUE} increases considerably the fitting time especially if surfaces are fitted i.e. \code{lo(~x1+x2)}.}
\item{obj}{an \code{lowss} object fitted within gamlss}
 \item{se}{if \code{se>0} then standard errors surfaces are drawn in the 3-dimensional plot. Set \code{se} at the required level i.e \code{se=1.96} will be an approximated 95\% CI.}
 \item{rug}{whether to plot a rug in the plot}
 \item{partial.resid}{whether to plot the partial residuals} 
 \item{col.term}{the colour of the line of fitted term}
 \item{cex.res}{ the shading of standard}
 \item{col.shaded}{the shading of standard error intervals}
 \item{col.res}{the colour of partial residuals}
 \item{col.rug}{the colour of the rug}
 \item{lwd.term}{the width of the line}
 \item{pch.res}{The character for the partial residuals}
 \item{type}{The type of the plot if the x's are two dimensional}
 \item{col.surface}{the colour of the fitted surface}
 \item{nlevels}{the number of levels used in \code{cotour()} plot.}
 \item{n.grid}{The number of points to evaluate the surface} 
 \item{image}{whether to use \code{image()} or just \code{contour}} 
}
\details{
 Note that \code{lo} itself does no smoothing; it simply sets things up for the function \code{gamlss.lo()} which is used by the backfitting function  \code{gamlss.add()}.
}
\value{
 a \code{loess} object is returned.  
}
\references{
Chambers, J. M. and Hastie, T. J. (1991). \emph{Statistical Models in S}, Chapman and Hall, London. 

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

\author{Mikis Stasinopoulos, Bob Rigby, (The original \code{lo()} function was based on the Trevor Hastie's S-plus \code{lo()} function. See also the documentation of the \code{loess} function for the authorship of the function.}
\note{Note that \code{lo} itself does no smoothing; it simply sets things up for \code{gamlss.lo()} to do the backfitting.}

\section{Warning }{In this version the first argument is a formula NOT a list as in the previous one} 

\seealso{ \code{\link{cs}}, \code{\link{random}},  }

\examples{
# fitting a loess curve with span=0.4 plus the a quarterly  effect 
aids1<-gamlss(y~lo(~x,span=0.4)+qrt,data=aids,family=PO) # 
term.plot(aids1, page=1)
\dontrun{
r1 <- gamlss(R~lo(~Fl)+lo(~A), data=rent, family=GA)
 term.plot(r1, pages=1)
 vis.lo(getSmo(r1, which=1), partial=T)
r2 <- gamlss(R~lo(~Fl+A), data=rent, family=GA)
 term.plot(r2, pages=1)
 vis.lo(getSmo(r2, which=1))
 vis.lo(getSmo(r2, which=1), se=1.97)
 vis.lo(getSmo(r2, which=1), partial.res=T)
}
}
\keyword{regression}% 
