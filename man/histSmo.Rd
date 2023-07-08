\name{histSmo}
\alias{histSmo}
\alias{histSmoC}
\alias{histSmoO}
\alias{histSmoP}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Density estimation using the Poisson trick
}
\description{
This set of functions use the old Poisson trick of discretising the data and then fitting a Poisson error model to the resulting frequencies (Lindsey, 1997). Here the model fitted is a smooth cubic spline curve. The result is a density estimator for the data.
}
\usage{
histSmo(y, lambda = NULL, df = NULL, order = 3, lower = NULL,  
       upper = NULL, type = c("freq", "prob"), 
       plot = FALSE, breaks = NULL,  
       discrete = FALSE, ...)
histSmoC(y, df = 10, lower = NULL, upper = NULL, type = c("freq", "prob"), 
       plot = FALSE, breaks = NULL,  
       discrete = FALSE, ...)
histSmoO(y, lambda = 1, order = 3, lower = NULL, upper = NULL, 
      type = c("freq", "prob"), 
      plot = FALSE, breaks = NULL,  
      discrete = FALSE, ...)
histSmoP(y, lambda = NULL, df = NULL, order = 3, lower = NULL, 
      upper = NULL, type = c("freq", "prob"), 
      plot = FALSE, breaks = NULL,  discrete = FALSE, 
      ...) 
}

\arguments{
  \item{y}{the variable of interest}
  \item{lambda}{the smoothing parameter}
  \item{df}{the degrees of freedom}
  \item{order}{the order of the P-spline}
  \item{lower}{the lower limit of the y-variable}
  \item{upper}{the upper limit of the y-variable}
  \item{type}{the type of histogram}
  \item{plot}{whether to plot the resulting density estimator}
  \item{breaks}{the number of break points to be used in the histogram and consequently the number of observations in the Poisson fit} 
  \item{discrete}{whether to treat the fitting density as a discrete distribution or not}
  \item{\dots}{further arguments passed to or from other methods.} 
}
\details{
Here are the methods used here:

i) The function \code{histSmoO()} uses Penalised discrete splines (Eilers,  2003). This function is appropriate when the smoothing parameter is fixed. 

ii) The function \code{histSmoC()} uses smooth cubic splines and fits a Poison error model to the frequencies using the \code{cs()} additive function of GAMLSS. This function is appropriate if the effective degrees of freedom are fixed in the model.

iii) The function \code{histSmoP()} uses Penalised cubic splines (Eilers and Marx 1996). It is fitting a Poisson model to the frequencies using the \code{pb()} additive function of GAMLSS. This function is appropriate if automatic selection of the smoothing parameter is required.

iv) The function \code{histSmo()} combines all the above functions in the sense that if lambda is fixed it uses  \code{histSmoO()}, if the df's are fixed it uses code{histSmoC()} and if none of these is specified it uses \code{histSmoP()}.

}
\value{
Returns a \code{histSmo} S3 object. The object has the following components:  
  \item{x}{the middle points of the discretise data}
  \item{counts}{how many observation are on the discretise intervals}
  \item{density}{the density value for each discrete interval}
  \item{hist}{the \code{hist} object used to discretise the data}
  \item{cdf}{The resulting cumulative distribution function useful for calculating probabilities from the estimate density}
  \item{nvcdf}{The inverse cumulative distribution function}
  \item{model}{The fitted Poisson model only for \code{histSmoP()} and \code{histSmoC()} }
}
\references{

Eilers, P. (2003). A perfect smoother. \emph{Analytical Chemistry}, 75: 3631-3636.

Eilers, P. H. C. and Marx, B. D. (1996). Flexible smoothing with
B-splines and penalties (with comments and rejoinder). \emph{Statist. Sci}, \bold{11}, 89-121.

Lindsey, J.K. (1997) \emph{Applying Generalized Linear Models}. New York: Springer-Verlag. 
ISBN 0-387-98218-3 

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
\author{
Mikis Stasinopoulos, Paul Eilers, Bob Rigby and Vlasios Voudouris
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{pb}}, \code{\link{cs}}
}
\examples{
# creating data from Pareto 2 distribution
set.seed(153)
Y <- rPARETO2(1000) 
\dontrun{
# getting the density 
histSmo(Y, lower=0, plot=TRUE)
# more  breaks a bit slower
histSmo(Y, breaks=200, lower=0, plot=TRUE)
# quick fit using lambda
histSmoO(Y, lambda=1, breaks=200, lower=0, plot=TRUE)
# or 
histSmo(Y, lambda=1, breaks=200, lower=0, plot=TRUE)
# quick fit using df
histSmoC(Y, df=15, breaks=200, lower=0,plot=TRUE)
# or 
histSmo(Y, df=15, breaks=200, lower=0)
# saving results
m1<- histSmo(Y, lower=0, plot=T)
plot(m1)
plot(m1, "cdf")
plot(m1, "invcdf")
# using with a histogram
library(MASS)
truehist(Y)
lines(m1, col="red")
#---------------------------
# now gererate from SHASH distribution
YY <- rSHASH(1000)
m1<- histSmo(YY)
# calculate Pr(YY>10)
1-m1$cdf(10)
# calculate Pr(-10<YY<10)
1-(1-m1$cdf(10))-m1$cdf(-10)
#---------------------------
#   from discrete distribution
YYY <- rNBI(1000, mu=5, sigma=4)
histSmo(YYY, discrete=TRUE, plot=T)
#
YYY <- rPO(1000, mu=5)
histSmo(YYY, discrete=TRUE, plot=T)
#
YYY <- rNBI(1000, mu=5, sigma=.1)
histSmo(YYY, discrete=TRUE, plot=T)
# genarating from beta distribution
YYY <- rBE(1000, mu=.1, sigma=.3)
histSmo(YYY, lower=0, upper=1, plot=T)
# from trucated data
Y <- with(stylo, rep(word,freq))
histSmo(Y, lower=1, discrete=TRUE, plot=T)
histSmo(Y, lower=1, discrete=TRUE, plot=T, type="prob")}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{distribution}
