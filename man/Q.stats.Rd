\name{Q.stats}
\alias{Q.stats}

\title{A function to calculate the Q-statistics}
\description{
 This function calculates and prints the Q-statistics (or Z-statistics) which are useful to test normality of the residuals within a range of an independent variable, for example age in centile estimation, see Royston and Wright (2000). 
}
\usage{
Q.stats(obj = NULL, xvar = NULL, resid = NULL, xcut.points = NULL, n.inter = 10, 
      zvals = TRUE, save = TRUE, plot = TRUE, digits.xvar = getOption("digits"),
      ...)
}

\arguments{
  \item{obj}{a GAMLSS object}
  \item{xvar}{a unique explanatory variable}
  \item{resid}{quantile or standardised residuals can be given here instead of a GAMLSS object in \code{obj}. In this case the function behaves differently (see details below) }
  \item{xcut.points}{the x-axis cut off points e.g. \code{c(20,30)}. If \code{xcut.points=NULL} then the \code{n.inter} argument is activated }
  \item{n.inter}{if \code{xcut.points=NULL} this argument gives the number of intervals in which the x-variable will be split, with default 10}
  \item{zvals}{if \code{TRUE} the output matrix contains the individual Z-statistics rather that  the Q statistics}
  \item{save}{whether to save the Q-statistics or not with default equal to \code{TRUE}. 
            In this case the functions produce a matrix giving individual Q (or z) statistics and the final aggregate Q's}
 \item{plot}{whether to plot  a visual version of the Q statistics  (default is TRUE)  }
  \item{digits.xvar}{to control the number of digits of the \code{xvar} in the plot }
   \item{\dots}{for extra arguments}
}
\details{
Note that the function \code{Q.stats} behaves differently depending whether the \code{obj} or the \code{resid} argument is set. The \code{obj} argument produces the Q-statistics (or Z-statistics) table appropriate for centile estimation (therefore it expect a reasonable large number of observations). The argument \code{resid} allows any model residuals, (not necessary GAMLSS), suitable standardised and is appropriate for any size of data. The resulting table contains only the  individuals  Z-statistics. 
}

\value{
  A table containing the  Q-statistics or Z-statistics. If \code{plot=TRUE} it produces also an graphical representation of the table. 
}
\references{

Rigby, R. A. and  Stasinopoulos D. M. (2005). Generalized additive models for location, scale and shape,(with discussion), 
\emph{Appl. Statist.}, \bold{54}, part 3, pp 507-554.

Royston P. and Wright E. M. (2000) Goodness of fit statistics for the age-specific reference intervals. 
\emph{Statistics in Medicine}, 19, pp 2943-2962.  

Rigby, R. A., Stasinopoulos, D. M.,  Heller, G. Z.,  and De Bastiani, F. (2019)
	\emph{Distributions for modeling location, scale, and shape: Using GAMLSS in R}, Chapman and Hall/CRC. An older version can be found in \url{https://www.gamlss.com/}.

Stasinopoulos D. M. Rigby R.A. (2007) Generalized additive models for location scale and shape (GAMLSS) in R.
\emph{Journal of Statistical Software}, Vol. \bold{23}, Issue 7, Dec 2007, \url{https://www.jstatsoft.org/v23/i07/}.

Stasinopoulos D. M., Rigby R.A., Heller G., Voudouris V., and De Bastiani F., (2017)
\emph{Flexible Regression and Smoothing: Using GAMLSS in R},  Chapman and Hall/CRC.  

(see also \url{https://www.gamlss.com/}).
}
\author{Mikis Stasinopoulos \email{d.stasinopoulos@londonmet.ac.uk}, Bob Rigby with contributions from Elaine Borghie}

\seealso{  \code{\link{gamlss}}, \code{\link{centiles.split}},  \code{\link{wp}} }
\examples{
data(abdom)
h<-gamlss(y~pb(x), sigma.formula=~pb(x), family=BCT, data=abdom) 
Q.stats(h,xvar=abdom$x,n.inter=8)
Q.stats(h,xvar=abdom$x,n.inter=8,zvals=FALSE)
Q.stats(resid=resid(h),  xvar=abdom$x, n.inter=5)
rm(h)
}
\keyword{regression}% 
