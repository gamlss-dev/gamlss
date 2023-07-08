\name{rvcov}
\alias{rvcov}
\alias{get.K}
\title{
Robust Variance-Covariance matrix of the parameters from a fitted GAMLSS model
}
\description{
The function \code{rvcov()} is design for providing  robust standard errors for the parameters estimates of a GAMLSS fitted model. The same result can be achieved by using \code{vcov(fitted_model,robust=TRUE)}. The function \code{get.()} gets the \eqn{K} matrix (see details below).  
}
\usage{
rvcov(object, type = c("vcov", "cor", "se", "coef", "all"),
                      hessian.fun = c("R", "PB") )
get.K(object, what = c("K", "Deriv"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
a GAMLSS fitted object
}
  \item{type}{this argument for \code{rvcov()} function 
whether variance-covariance matrix, correlation matrix, standard errors  or  all of them  
}
  \item{what}{this an argument for the function \code{ket.K()} allowing to get either \eqn{K} or the first derivative of the likelihood with respect to the parameters (the \eqn{\beta}'s in  the GAMLSS notation).} 
\item{hessian.fun}{How to obtain numerically the Hessian i) using \code{optimHess()}, option     \code{"R"}  ii) using a function by Pinheiro and Bates taken from package \code{nlme}, option \code{"PB"}.}
}

\details{
The robust standard errors are calculated for the robust sandwich estimator of the variance-covariance given by \eqn{S=VKV} where \eqn{V} is the standard variance-covariance matrix (the inverse of the information matrix) and \eqn{K} is an estimate of the variance of he first derivatives of he likelihood.  The function \code{get.K()} is use the get the required \eqn{K} matrix.
}
\value{
A variance covariance matrix or other relevant output
}
\references{
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
\author{Mikis Stasinopoulos \email{d.stasinopoulos@londonmet.ac.uk}, Bob Rigby and Vlasios Voudouris}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{vcov}}, ~~~
}
\examples{
# gererate from a gamma distribution 
Y <- rGA(200, mu=1, sigma=2)
hist(Y)
# fitting the wrong model i.e. sigma=1 
m1 <- gamlss(Y~1, family=EXP) 
# the conventinal se is too precise 
vcov(m1, type="se")
# the sandwich se is wider  
rvcov(m1, type="se") 
# fitting the correct model 
 m2 <- gamlss(Y~1, family=GA)
 vcov(m2, type="se")
 rvcov(m2, type="se")
# similar stadard errors
# also obtained using 
vcov(m2, type="se", robust=TRUE)


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{regression}% 

