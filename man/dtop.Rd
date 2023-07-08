\name{dtop}
\alias{dtop}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Detrended transformed Owen's plot}
\description{
Provides single or multiple detrended transformed  Owen's plot, Owen (1995), for a GAMLSS 
fitted objects or any other fitted object which has the method resid(). This is a diagnostic tool for checking whether the normalised 
quantile residuals are coming from a normal distribution or not. 
This could be true if the horizontal line is within the confidence intervals.
}
\usage{
dtop(object = NULL, xvar = NULL, resid = NULL,
      type = c("Owen", "JW"), 
      conf.level = c("95", "99"), n.inter = 4, 
      xcut.points = NULL, overlap = 0, 
      show.given = TRUE, cex = 1, pch = 21, 
      line = TRUE, ...)}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{a GAMLSS fitted object or any other fitted object which has the method resid().}
  \item{xvar}{the explanatory variable against which the detrended Owen's plots will be plotted}
  \item{resid}{if the object is not specified the residual vector can be given here}
  \item{type}{whether to use Owen (1995) or Jager and Wellner (2004) approximate formula}
  \item{conf.level}{95 (default) or 99 percent confidence interval for the plots}
  \item{n.inter}{he number of intervals in which the explanatory variable xvar will be cut}
  \item{xcut.points}{the x-axis cut off points e.g. c(20,30). If xcut.points=NULL then the n.inter argument is activated}
  \item{overlap}{how much overlapping in the xvar intervals. Default value is overlap=0 for non overlapping intervals}
  \item{show.given}{whether to show the x-variable intervals in the top of the graph, default is show.given=TRUE}
  \item{cex}{the cex plotting parameter with default cex=1}
  \item{pch}{the pch plotting parameter with default pch=21}
  \item{line}{whether the detrended empirical cdf should be plotted or not}
  \item{\dots}{for extra arguments}
}
\details{If the xvar argument is not specified then a single detrended Owen's plot is used, see Owen (1995). In this case the plot is a detrended nonparametric likelihood confidence band for a distribution function. That is, if the horizontal lines lies within the confidence band then the normalised residuals could have come from a Normal distribution and consequently the assumed response variable distribution is reasonable. If the xvar is specified then we have as many plots as n.iter. In this case the x-variable is cut into n.iter intervals with an equal number observations and detrended Owen's plots for each interval are plotted. This is a way of highlighting failures of the model within different ranges of the explanatory variable.}
\value{A plot is returned.}
\references{

Jager, L. and Wellner, J. A (2004) A new goodness of fit test: the reversed Berk-Jones statistic, University of Washington, Department of Statistics, Technical report 443.    

Owen A. B. (1995) Nonparametric Confidence Bands for a Distribution Function. Journal of the American Statistical Association Vol. 90, No 430, pp. 516-521.

Rigby, R. A. and Stasinopoulos D. M. (2005). Generalized additive models for location, scale and shape,(with discussion), Appl. Statist., 54, part 3, 1-38.

Rigby, R. A., Stasinopoulos, D. M.,  Heller, G. Z.,  and De Bastiani, F. (2019)
	\emph{Distributions for modeling location, scale, and shape: Using GAMLSS in R}, Chapman and Hall/CRC. An older version can be found in \url{https://www.gamlss.com/}.

Stasinopoulos D. M. Rigby R.A. (2007) Generalized additive models for location scale and shape (GAMLSS) in R.
\emph{Journal of Statistical Software}, Vol. \bold{23}, Issue 7, Dec 2007, \url{https://www.jstatsoft.org/v23/i07/}.

Stasinopoulos D. M., Rigby R.A., Heller G., Voudouris V., and De Bastiani F., (2017)
\emph{Flexible Regression and Smoothing: Using GAMLSS in R},  Chapman and Hall/CRC.  

(see also \url{https://www.gamlss.com/}).
}
\author{Mikis Stasinopoulos, Bob Rigby and Vlassios Voudouris}
\seealso{\code{\link{wp}}}
\examples{
data(abdom)
a<-gamlss(y~pb(x),sigma.fo=~pb(x,1),family=LO,data=abdom)
dtop(a)
dtop(a, xvar=abdom$x)
rm(a)
}
\keyword{regression}% 

