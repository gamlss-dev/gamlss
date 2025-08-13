\name{acfResid}
\alias{acfResid}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
ACF plot of the residuals
}
\description{
This plot display the ACF and PACF of the residuals of a gamlss or other fitted model (provided that they have been standardised appropriately.  Is is appropriate for time series data.
}
\usage{
acfResid(obj = NULL, resid = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obj}{A gamlss model or other fitted model where the resid() function applies exist
}
  \item{resid}{
if \code{obj} does not exist the argument here will be used
}
}
\details{
The ACF and PACF for the residuals \code{r}, squared residuals \code{r^2},  \code{r^3} and \code{r^4} are plotted
}
\value{
The relevant plots are displayed
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


}
\author{
Mikis Stasinopoulos. Bob Rigby. Vlasios Voudouris and Majid Djennad
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{acf}}
}
\examples{
library(datasets)
data(co2)
m1<- gamlss(co2~pb(as.numeric(time(co2)))+factor(as.numetic(cycle(co2))))
acfResid(m1)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{regresson}
\keyword{ts}% __ONLY ONE__ keyword per line
