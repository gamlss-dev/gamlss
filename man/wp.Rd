\name{wp}
\alias{wp}
\title{ Worm plot}
\description{
Provides a single plot or multiple worm plots for a GAMLSS fitted or more general for any fitted models where the method \code{resid()} exist and the residuals are defined sensibly. The worm plot (a de-trended QQ-plot), van Buuren and Fredriks M. (2001), is a diagnostic tool for checking the residuals within different ranges (by default not overlapping) of the explanatory variable(s).
          }
\usage{
wp(object = NULL, xvar = NULL, resid = NULL, n.inter = 4,
  xcut.points = NULL, overlap = 0, xlim.all = 4, 
  xlim.worm = 3.5, show.given = TRUE, line = TRUE, 
  ylim.all = 12 * sqrt(1/length(resid)), 
  ylim.worm = 12 * sqrt(n.inter/length(resid)), 
  cex = 1, cex.lab = 1, pch = 21, ...)
}

\arguments{
  \item{object}{a GAMLSS fitted object or any other fitted model where the \code{resid()} method works (preferably it should be standardised or quantile residuals)}
  \item{xvar}{the explanatory variable(s) against which the worm plots will be plotted. If only one variable is involved use \code{xvar=x1} if two variables are involved use \code{xvar=~x1*x2}. See also note below for use of formula if the data argument is not found in the fitted model}
  \item{resid}{if object is missing this argument can be used to specify the residual vector (again it should a quantile residuals or it be assumed to come from a normal distribution)}
  \item{n.inter}{the number of intervals in which the explanatory variable \code{xvar} will be cut}
  \item{xcut.points}{the x-axis cut off points e.g. \code{c(20,30)}. If \code{xcut.points=NULL} then the \code{n.inter} argument is activated }
  \item{overlap}{how much overlapping in the \code{xvar} intervals. Default value is \code{overlap=0} for non overlapping intervals}
  \item{xlim.all}{for the single plot, this value is the x-variable limit, default is \code{xlim.all=4}}
  \item{xlim.worm}{for multiple plots, this value is the x-variable limit, default is \code{xlim.worm=3.5}}
  \item{show.given}{whether to show the x-variable intervals in the top of the graph, default is \code{show.given=TRUE} }
  \item{line}{whether to plot the polynomial line in the worm plot, default value is \code{line=TRUE}}
  \item{ylim.all}{for the single plot, this value is the y-variable limit, default value is \code{ylim.all=12*sqrt(1/length(fitted(object)))}}
  \item{ylim.worm}{for multiple plots, this values is the y-variable limit, default value is \code{ylim.worm=12*sqrt(n.inter/length(fitted(object)))}}
  \item{cex}{ the cex plotting parameter for changing the side of worm with default \code{cex=1}}
  \item{cex.lab}{the cex plotting parameter for changing the size of the axis labels}
  \item{pch}{ the pch plotting parameter with default \code{pch=21} }
   \item{\dots}{for extra arguments}
 }
\details{
 If the \code{xvar} argument is not specified then a single worm plot is used. In this case a worm plot is a de-trended normal QQ-plot so departure from normality is highlighted. 

If a single  \code{xvar} is specified (with or without the use of a formula) i.e. \code{xvar=x1} or \code{xvar=~x1}) then we have as many worm plot as \code{n.iter}. 
 In this case the x-variable is cut into \code{n.iter} intervals with an equal number observations and de-trended normal QQ (i.e. worm) plots for each interval are plotted. This is a way of highlighting failures of the model within different ranges of the 
the single explanatory variable. The fitted coefficients from fitting cubic polynomials to the residuals (within each x-variable interval) can be obtain by e.g. \code{coeffs<-wp(model1,xvar=x,n.iner=9)}.  van Buuren and Fredriks M. (2001) used these residuals to identify regions (intervals) of the explanatory variable within which the model does not fit adequately the data (called "model violation")  

Two variables can be displayed with the use of a formula, i.e. \code{xvar=~x1*x2}. In this case the \code{n.inter} can be a vector with two values.   
 
}
\value{
  For multiple plots the \code{xvar} intervals and the coefficients of the fitted cubic polynomials to the residuals (within each \code{xvar} interval) are returned.   
}
\note{Note that the \code{wp()} function, if the argument \code{object} is used, is looking for the data argument of the object. If the argument \code{data} exists it uses its environment to find   \code{xvar} (whether it is a formula or not). As a result if \code{data} exists within \code{object}  
\code{xvar=~x*f} can be used (assuming that \code{x} and \code{f} are in the data) otherwise the variable should be explicitly defined i.e. \code{xvar=~data$x*data$f}.
}
\references{
Rigby, R. A. and  Stasinopoulos D. M. (2005). Generalized additive models for location, scale and shape,(with discussion), \emph{Appl. Statist.}, \bold{54}, part 3,
1-38. 
 
Stasinopoulos D. M. Rigby R.A. (2007) Generalized additive models for location scale and shape (GAMLSS) in R.
\emph{Journal of Statistical Software}, Vol. \bold{23}, Issue 7, Dec 2007, \url{http://www.jstatsoft.org/v23/i07}.

Stasinopoulos D. M., Rigby R.A., Heller G., Voudouris V., and De Bastiani F., (2017) \emph{Flexible Regression and Smoothing: Using GAMLSS in R}, Chapman and Hall/CRC. 
(see also \url{http://www.gamlss.org/}).
           
van Buuren and Fredriks M. (2001) Worm plot: simple diagnostic device for modelling growth reference curves. 
            \emph{Statistics in Medicine}, \bold{20}, 1259--1277
            }
\author{Mikis Stasinopoulos and Bob Rigby}


\seealso{  \code{\link{gamlss}}, \code{\link{plot.gamlss} }}

\examples{
data(abdom)
# with data
a<-gamlss(y~pb(x),sigma.fo=~pb(x,1),family=LO,data=abdom)
wp(a)
coeff1<-wp(a,xvar=x)
coeff1
\dontrun{
# no data argument
b <- gamlss(abdom$y~pb(abdom$x),sigma.fo=~pb(abdom$x),family=LO)
wp(b) 
wp(b, xvar=abdom$x)# not wp(b, xvar=x)
# using  the argument resid
# this will work
wp(resid=resid(a),  xvar=abdom$x)
# not this
# wp(resid=resid(a),  xvar=x)
# this example uses the rent data
m1 <- gamlss(R~pb(Fl)+pb(A)+loc, sigma.fo=~pb(Fl)+pb(A), data=rent, family=GA)
# a single worm plot
wp(m1, ylim.all=0.5)
# a single continuous x variable 
wp(m1, xvar=Fl, ylim.worm=.8)
# a single x variable changing the default number of intervals
wp(m1, xvar=Fl, ylim.worm=1.5, n.inter=9)
# different x variable changing the default number of intervals
B1<-wp(m1, xvar=A, ylim.worm=1.2, n.inter=9) 
B1
# the number five plot has intervals 
# [5,] 1957.5 1957.5 
# rather disappoining 
# try formula for xvar
wp(m1, xvar=~A, ylim.worm=1.2, n.inter=9)
# better in this case using formula
# now using a factor included in the model
wp(m1, xvar=~loc, ylim.worm=1.2, n.inter=9)
# using a factor notin the model
wp(m1, xvar=~B, ylim.worm=1.5, n.inter=9)
# level 2 (with B=1) did not fit well
# trying two continuous variable 
wp(m1, xvar=~Fl*A, ylim.worm=1.5, n.inter=4)
# one continuous and one categorical 
wp(m1, xvar=~Fl*loc, ylim.worm=1.5, n.inter=4)
# two categorical
wp(m1, xvar=~B*loc, ylim.worm=1.5, n.inter=4)
}

}
\keyword{regression}% 
