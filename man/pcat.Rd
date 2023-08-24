\name{pcat}
\alias{pcat}
\alias{gamlss.pcat}
\alias{plotLambda}
\alias{plotDF}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Reduction for the Levels of a Factor.
}
\description{
The function is trying to merged similar levels of a given factor. Its based on ideas given by 
Tutz (2013). 
}
\usage{
pcat(fac, df = NULL, lambda = NULL, method = c("ML", "GAIC"), start = 0.001, 
         Lp = 0, kappa = 1e-05, iter = 100, c.crit = 1e-04, k = 2)

gamlss.pcat(x, y, w, xeval = NULL, ...)

plotDF(y, factor = NULL, formula = NULL, data, along = seq(0, nlevels(factor)), 
         kappa = 1e-06, Lp = 0, ...)

plotLambda(y, factor = NULL, formula = NULL, data, along = seq(-2, 2, 0.1), 
         kappa = 1e-06, Lp = 0, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fac, factor}{a factor to reduce its levels} 
  \item{df}{the effective degrees of freedom df}
  \item{lambda}{the smoothing parameter}
  \item{method}{which method is used for the estimation of the smoothing parameter, \code{"ML"} or \code{"GAIC"} are allowed.}
  \item{start}{ starting value for \code{lambda} if it estimated using \code{"ML"} or \code{"GAIC"}}
  \item{Lp}{
The type of penalty required, \code{Lp=0} is the default. Use \code{Lp=1} for lasso type and different values for different required penalty.
}
  \item{kappa}{a regulation parameters used for the weights in the penalties.}
  \item{iter}{
the number of internal iteration allowed
}
  \item{c.crit}{the convergent criterion}
  \item{k}{the penalty if \code{"GAIC"} method is used.}
  \item{x}{explanatory factor}
  \item{y}{the response or iterative response variable}
  \item{w}{iterative weights}
  \item{xeval}{indicator whether to predict}
   \item{formula}{A formula}
    \item{data}{A data frame}
    \item{along}{a sequence of values} 
  \item{\dots}{for extra variables}
}
\details{
The \code{pcat()} is used for the fitting of the factor. The function shrinks the levels of the categorical factor (not towards the overall mean as the function \code{random()} is doing) but towards each other.  This results to a reduction of the number if levels of the factors.  Different norms can be used for the shrinkage by specifying the argument \code{Lp}.  
}
\value{
The function \code{pcat} reruns a vector  endowed with a number of attributes. 
The vector itself is used in the construction of the model matrix, while the attributes are needed for the backfitting algorithms additive.fit(). The backfitting is done in \code{gamlss.pcat}.
}


\references{
Tutz G. (2013) Regularization and Sparsity in Discrete Structures in the \emph{Proceedings of the 29th International Workshop on Statistical Modelling}, Volume 1, p 29-42, Gottingen, Germany

Rigby, R. A., Stasinopoulos, D. M.,  Heller, G. Z.,  and De Bastiani, F. (2019)
	\emph{Distributions for modeling location, scale, and shape: Using GAMLSS in R}, Chapman and Hall/CRC. An older version can be found in \url{https://www.gamlss.com/}.

Stasinopoulos D. M. Rigby R.A. (2007) Generalized additive models for location scale and shape (GAMLSS) in R.
\emph{Journal of Statistical Software}, Vol. \bold{23}, Issue 7, Dec 2007, \url{https://www.jstatsoft.org/v23/i07/}.

Stasinopoulos D. M., Rigby R.A., Heller G., Voudouris V., and De Bastiani F., (2017)
\emph{Flexible Regression and Smoothing: Using GAMLSS in R},  Chapman and Hall/CRC.  

(see also \url{https://www.gamlss.com/}).
}
\author{Mikis Stasinopoulos,  Paul Eilers and Marco Enea}
\note{Note that \code{pcat} itself does no smoothing; it simply sets things up for \code{gamlss.pcat()} to do the smoothing within the backfitting.}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{ \code{\link{random}} }
\examples{
# Simulate data 1
    n <- 10  # number of levels 
    m <- 200 # number of observations  
set.seed(2016)
level <-  as.factor(floor(runif(m) * n) + 1)
  a0  <-  rnorm(n)
sigma <-  0.4
   mu <-  a0[level]
   y <-  mu + sigma * rnorm(m)
plot(y~level)
points(1:10,a0, col="red")
 da1 <- data.frame(y, level)
#------------------
  mn <- gamlss(y~1,data=da1 ) # null model 
  ms <- gamlss(y~level-1, data=da1) # saturated model 
  m1 <- gamlss(y~pcat(level), data=da1) # calculating lambda ML
AIC(mn, ms, m1)
\dontrun{
m11 <- gamlss(y~pcat(level, method="GAIC", k=log(200)), data=da1) # GAIC
AIC(mn, ms, m1, m11) 
#gettng the fitted object -----------------------------------------------------
getSmo(m1)
coef(getSmo(m1))
fitted(getSmo(m1))[1:10]
plot(getSmo(m1)) # 
# After the fit a new factor is created  this factor has the reduced levels
 levels(getSmo(m1)$factor)
# -----------------------------------------------------------------------------
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{regeression}

