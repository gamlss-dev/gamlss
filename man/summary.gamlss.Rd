\name{summary.gamlss}
\alias{summary.gamlss}

\title{Summarizes a GAMLSS fitted model}
\description{\code{summary.gamlss} is the GAMLSS specific method for the generic function \code{summary} which summarize 
      objects returned by modelling functions. 
}
\usage{
\method{summary}{gamlss}(object, type = c("vcov", "qr"), 
                       robust=FALSE, save = FALSE,
                       hessian.fun = c("R", "PB"), 
                       digits = max(3, getOption("digits") - 3),...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{a GAMLSS fitted model}
  \item{type}{the default value \code{vcov} uses the \code{vcov()} method for gamlss to get the 
  variance-covariance  matrix of the estimated beta coefficients, see details below. 
  The alternative \code{qr} is the original method used in gamlss to 
  estimated the standard errors but it is not reliable since it do not take into the account the inter-correlation between
  the distributional parameters \code{mu}, \code{sigma}, \code{nu} and \code{tau}.
 }
  \item{robust}{whether robust (sandwich) standard errors are required}
  \item{save}{whether to save the environment of the function so to have access to its values}
  \item{hessian.fun}{whether when calculate the Hessian should use the "R" function \code{optimHess()} or a function based on Pinheiro and Bates \code{nlme} package, "PB".}
 \item{digits}{the number of digits in the output}
\item{\dots}{for extra arguments}
}
\details{Using the  default value \code{type="vcov"}, the \code{vcov()} method for gamlss is used to get  
  the variance covariance matrix (and consequently the standard errors) of the beta parameters. 
  The variance covariance matrix is  calculated using the inverse of the numerical second derivatives
  of the observed information matrix. This is a more reliable method since it take into the account the 
  inter-correlation between the all the parameters. The \code{type="qr"} assumes that the parameters are fixed 
  at the estimated values. Note that both methods are not appropriate and should be used with caution if smoothing 
  terms are used in the fitting. }
 
\value{
 Print summary of a GAMLSS object
}
\references{Rigby, R. A. and  Stasinopoulos D. M. (2005). Generalized additive models for location, scale and shape,(with discussion), 
\emph{Appl. Statist.}, \bold{54}, part 3, pp 507-554.


Rigby, R. A., Stasinopoulos, D. M.,  Heller, G. Z.,  and De Bastiani, F. (2019)
	\emph{Distributions for modeling location, scale, and shape: Using GAMLSS in R}, Chapman and Hall/CRC. An older version can be found in \url{https://www.gamlss.com/}.

Stasinopoulos D. M. Rigby R.A. (2007) Generalized additive models for location scale and shape (GAMLSS) in R.
\emph{Journal of Statistical Software}, Vol. \bold{23}, Issue 7, Dec 2007, \url{https://www.jstatsoft.org/v23/i07/}.

Stasinopoulos D. M., Rigby R.A., Heller G., Voudouris V., and De Bastiani F., (2017)
\emph{Flexible Regression and Smoothing: Using GAMLSS in R},  Chapman and Hall/CRC.  

(see also \url{https://www.gamlss.com/}).
 }
\author{Mikis Stasinopoulos, Bob Rigby  and Calliope Akantziliotou }

\seealso{  \code{\link{gamlss}}, \code{\link{deviance.gamlss}},  \code{\link{fitted.gamlss}} }
\examples{
data(aids)
h<-gamlss(y~poly(x,3)+qrt, family=PO, data=aids) # 
summary(h)
rm(h)
}
\keyword{regression}% 
