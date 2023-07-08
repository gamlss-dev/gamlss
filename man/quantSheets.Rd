\name{quantSheets}
\alias{quantSheets}
\alias{quantSheets.control}
\alias{findPower}
\alias{z.scoresQS}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Quantile Sheets
}
\description{
The quantile sheets function \code{quantSheets()} is based on the work of Sabine
Schnabe and Paul  Eiler (see references below). The estimation  of the quantile curves 
is done simultaneously by also smoothing in the direction of y as well as x. This avoids (but do not eliminate completely) the problem of crossing quantiles.
}
\usage{
quantSheets(y, x, x.lambda = 1, p.lambda = 1, data = NULL, 
            cent = 100 * pnorm((-4:4) * 2/3), 
            control = quantSheets.control(...), print = TRUE,  ...)

quantSheets.control(x.inter = 10, p.inter = 10, degree = 3, logit = FALSE, 
            order = 2, kappa = 0, n.cyc = 100, c.crit = 1e-05, plot = TRUE, 
            power = NULL, ...)

findPower(y, x, data = NULL, lim.trans = c(0, 1.5), prof = FALSE, 
            k = 2, c.crit = 0.01, step = 0.1)

z.scoresQS(object, y, x, plot = FALSE, tol = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{the y variable}
  \item{x}{the x variable}
  \item{x.lambda}{smoothing parameter in the direction of x}
  \item{p.lambda}{smoothing parameter in the direction of y (probabilities)}
  \item{data}{the data frame}
  \item{cent}{the centile values where the quantile sheets is evaluated}
  \item{control}{for the parameters controlling the algorithm}
  \item{print}{whether to print the sample percentages}
    \item{x.inter}{number of intervals in the x direction for the B-splines}
  \item{p.inter}{number of intervals in the probabilities (y-direction) for the B-splines}
  \item{degree}{the degree for the B-splines}
  \item{logit}{whether to use \code{logit(p)} instead of \code{p} (probabilities) for the y-axis}
  \item{order}{the order of the penalty}
  \item{kappa}{is a ridge parameter set to zero (for no ridge effect) }
  \item{n.cyc}{number of cycles of the algorithm}
  \item{c.crit}{convergence criterion of the algorithm}
  \item{plot}{whether to plot the resulting quantile sheets}
  \item{power}{The value of the power transformation in the x axis if needed}
  \item{lim.trans}{the limits for looking for the power transformation 
          parameter using \code{findPower()} }
  \item{prof}{whether to use the profile GAIC or \code{optim()} to the parameter
            the power transformation}
  \item{k}{the  GAIC penalty}
  \item{step}{the steps for the profile GAIC if the argument \code{prof} of 
              \code{findPower()} is TRUE} 
  \item{object}{a fitted \code{quantSheets} object}
  \item{tol}{how far out from the range of the y variable should go for 
             estimating the distribution of y using the \code{flexDist()} 
             function}
  \item{\dots}{for further arguments}
}
\details{
The advantage  of quantile sheets is that they estimates simultaneously all the quantiles.
This almost eliminates the problem of crossing quantiles. The method is very fast and
useful for  exploratory  tool. The function needs two smoothing parameters. 
Those two parameters have to specified by the user. They are \emph{not} estimated automatically. They can be selected by visual inspection.


The disadvantages of quantile sheets comes from the fact that like all non-parametric techniques do not have a goodness of fit measure to change how good is the models and the residuals based diagnostics are not existence since it is difficult to define residuals in this set up.  

In this implementation we do provide residuals by using the  \code{flexDist()} 
function from package \pkg{gamlss.dist}. This is based on the idea that by 
knowing the quantiles of the distribution we can reconstruct non parametrically 
the distribution itself and this is what \code{flexDist()} is doing.
As a word of caution, such a construct is based on several assumptions and depends on 
several smoothing parameters.   Treat those residuals with caution. 
The same caution should apply to the function \code{z.scoresQS()}. 
}
\value{Using the function \code{quantSheets()} a \code{quantSheets} object is returned having the following  methods:
\code{print()}, \code{fitted()}, \code{predict()} and \code{resid()}. 

Using \code{findPower()} a single values of the power parameter is returned.

Using \code{z.scoresQS} a vector of z-scores is returned.

}
\references{
Schnabel, S.K. (2011) \emph{Expectile smoothing: new perspectives on asymmetric least squares. An application to life expectancy}, Utrecht University.

Schnabel, S. K and Eilers, P. H. C.(2013) Simultaneous estimation of quantile curves using quantile sheets, \emph{AStA Advances in Statistical Analysis}, \bold{97}, 1, 
  pp 77-87, Springer.


Schnabel, S. K and Eilers, P. H. (2013) A location-scale model for non-crossing 
expectile curves, \emph{Stat}, \bold{2}, 1, pp 171-183.

Rigby, R. A., Stasinopoulos, D. M.,  Heller, G. Z.,  and De Bastiani, F. (2019)
	\emph{Distributions for modeling location, scale, and shape: Using GAMLSS in R}, Chapman and Hall/CRC. An older version can be found in \url{https://www.gamlss.com/}.

Stasinopoulos D. M. Rigby R.A. (2007) Generalized additive models for location scale and shape (GAMLSS) in R.
\emph{Journal of Statistical Software}, Vol. \bold{23}, Issue 7, Dec 2007, \url{https://www.jstatsoft.org/v23/i07/}.

Stasinopoulos D. M., Rigby R.A., Heller G., Voudouris V., and De Bastiani F., (2017)
\emph{Flexible Regression and Smoothing: Using GAMLSS in R},  Chapman and Hall/CRC.  

(see also \url{https://www.gamlss.com/}).
 }

\author{Mikis Stasinopoulos based on function provided by Paul Eiler and Sabine
Schnabe}



\seealso{
\code{\link{lms}}: for  a parametric equivalent results.
}
\examples{
data(abdom)
m1 <- quantSheets(y,x, data=abdom)
head(fitted(m1))
p1 <- predict(m1, newdata=c(20,30,40))
matpoints(c(20,30,40), p1)
z.scoresQS(m1,y=c(150, 300),x=c(20, 30) )
# If we needed a power transformation not appropriate for this data
findPower(y,x, data=abdom)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{regression}
