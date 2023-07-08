\name{getQuantile}
\alias{getQuantile}

\title{Getting the partial quantile function for a term}
\description{
This function can be used to calculate the partial effect that an explanatory variable has on a specific quantile.

By `partial effect' function we mean how the term influence the quantile given that the rest of explanatory terms are constant.  

The function takes a GAMLSS object and for the range of a specified explanatory
(by fixing the rest of the terms at specified values), 
calculates the effect that this term has on the a quantile of the distribution.
That is, it shows the effect that the particular term  has on the quantile. 
The `partial' quantile is calculated on a finite grid of values and then the function is approximated (using the \code{splinefun()})  and  saved.

The saved function can be used to calculate the first derivative. This first derivatives   shows the chance of the quantile function  for a small change in the explanatory variable, by fixing the rest of the explanatory variables at a constant values.
}
\usage{
getQuantile(obj = NULL, term = NULL, quantile = 0.9, data = NULL,
    n.points = 100, how = c("median", "last"), 
    fixed.at = list(), plot = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obj}{A \code{gamlss} object}
  \item{term}{
an explanatory variable (at the moment works with with continuous)
}
  \item{quantile}{the required quantile of the distribution}
  \item{data}{the data.frame (not needed if is declared on \code{obj})}
  \item{n.points}{the number of points in which the quantile function needs evaluation}
  \item{how}{whether for extra continuous explanatory variables should fixed at the median or the last observation in the data}
  \item{fixed.at}{a list indicating at which values the rest of the explanatory terms should be fixed}
  \item{plot}{whether to the plot the partial quantile function and its first derivatives}
}
\details{
The function \code{getQuantile()} relies on the \code{predictAll()} function to evaluate the distribution parameters at a grid (default 100 points) of the specified term (given that the  the rest of the terms are fixed). Then the inverse cdf is used to calculate the partial quantile. The function then is approximated using   \code{splinefun()})  and saved. 
}
\value{
A function is created which can be used to evaluate the partial effect of the explanatory variable  on a specified quantile. 
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
\author{Mikis Stasinopoulos}

\seealso{\code{\link{gamlss}}, \code{\link{getPEF}}}
\examples{
library(gamlss)
data(rent)
m1 <- gamlss(R~pb(Fl)+pb(A)+B+loc, data=rent, family=GA)
FF<-getQuantile(m1, quantile=0.9, term="A", plot=TRUE)
FF(1960)
FF(1060, deriv=1)
FF(1060, deriv=2)
\dontrun{
# plotting partial quantile
# .05, 0.25, 0.5, 0.75, 0.95
# at the default values
# Fl = median(Fl), B=0, and loc=2
plot(R~A, data=rent,  col="lightgray", pch=20)
for (i in c(.05, 0.25, 0.5, 0.75, 0.95))
{
  Qua <- getQuantile(m1, quantile=i,term="A")
  curve(Qua, 1900, 1985,  lwd=1, lty=1, add=T)
}
# plotting at values Fl=60, B=1, and loc=1.
for (i in c(.05, 0.25, 0.5, 0.75, 0.95))
{
  Qua <- getQuantile(m1, quantile=i,term="A", 
            fixed.at=list(Fl=60, B=1, loc=1))
  curve(Qua, 1900, 1985,  lwd=1, lty=2, col="red", add=T)
}
# plotting at Fl=60, B=1 and loc=1.
for (i in c(.05, 0.25, 0.5, 0.75, 0.95))
{
  Qua <- getQuantile(m1, quantile=i,term="A", 
           fixed.at=list(Fl=120, B=0, loc=3))
  curve(Qua, 1900, 1985,  lwd=1, lty=3, col="blue", add=T)
}

}
}

\keyword{Statistical Models}
\keyword{regression}
