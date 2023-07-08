\name{calibration}
\alias{calibration}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Calibrating centile curves}
\description{
 This function can used when the fitted model centiles do not coincide with the sample centiles.  
}
\usage{
calibration(object, xvar, cent = c(0.4, 2, 10, 25, 50, 75, 90, 98, 99.6),
            legend = FALSE, fan = FALSE, ...)
}
\arguments{
  \item{object}{a gamlss fitted object}
  \item{xvar}{The explanatory variable}
  \item{cent}{a vector with elements the \% centile values for which the centile curves have to be evaluated}
  \item{legend}{whether legend is required}
  \item{fan}{whether to use the fan  version of centiles}
  \item{\dots}{other argument pass on to \code{centiles()} function}
}
\details{
The  function finds the sample quantiles of the residuals of the fitted model (the z-scores) and use them as sample quantile in the argument \code{cent} of the \code{centiles()} function. This procedure is appropriate if the fitted model centiles do not coincide with the sample centiles and when this failure is the same in all values of the explanatory variable \code{xvar}.
}
\value{
A centile plot is produced and the sample centiles below each centile curve are printed (or saved)
}
\references{
Rigby, R. A. and  Stasinopoulos D. M. (2005). Generalized additive models for location, scale and shape, (with discussion), 
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
Mikis Stasinopoulos \email{d.stasinopoulos@londonmet.ac.uk}, Bob Rigby \email{r.rigby@londonmet.ac.uk} and Vlasios Voudouris}

\seealso{\code{\link{centiles}}, \code{\link{centiles.fan}}}
\examples{
data(abdom)
m1<-gamlss(y~pb(x), sigma.fo=~pb(x), family=LO, data=abdom)
calibration(m1, xvar=abdom$x, fan=TRUE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{regression}% 
