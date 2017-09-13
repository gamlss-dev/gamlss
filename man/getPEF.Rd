\name{getPEF}
\alias{getPEF}

\title{Getting the partial effect function from a continuous term in a GAMLSS model}
\description{
This function can be used to calculate the partial effect and the elasticity of a continuous explanatory variable \code{x}.

By `partial effect' function we mean how \code{x} is influence the parameter of interest given that the rest of explanatory terms for this paramerer are on (specified) fixed values.  

The function takes a GAMLSS object and for the range of the continuous variable \code{x}, 
(by fixing the rest of the explanatory terms at specified values), 
calculates the effect that \code{x} has on the specific distribution parameter (or its predictor).
The resulting function shows the effect that \code{x} has on the distribution parameter. 
The partial effect function which is calculated on a finite grit is then approximated using the \code{splinefun()} in R and its is saved.

The saved function can be used to calculate the elasticity of \code{x}. The elasticity is the first derivative of the partial effect function  and shows the chance of the parameter of interest for a small change in in \code{x}, by fixing the rest of the explanatory variables at specified values.
}
\usage{
getPEF(obj = NULL, term = NULL, data = NULL, n.points = 100, 
          parameter = c("mu", "sigma", "nu", "tau"), 
          type = c("response", "link"), how = c("median", "last"), 
          fixed.at = list(), plot = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obj}{A \code{gamlss} object}
  \item{term}{the continuous explanatory variable}
  \item{data}{the data.frame (not needed if is declared on \code{obj})}
  \item{n.points}{the number of points in which the influence function for \code{x} need to be evaluated}
  \item{parameter}{which distribution parameter}
  \item{type}{whether against the parameter, \code{"response"}, or the predictor \code{"link"}}
  \item{how}{whether for continuous variables should use the median  or the last observation in the data}
  \item{fixed.at}{a list indicating  at which values the rest of the explanatory terms should be fixed}
  \item{plot}{whether to the plot the influence function and its first derivatives}
}

\value{A function is created which can be used to evaluate the partial effect function at different values of \code{x}.  }
\references{
Rigby, R. A. and  Stasinopoulos D. M. (2005). Generalized additive models for location, scale and shape,(with discussion), 
\emph{Appl. Statist.}, \bold{54}, part 3, pp 507-554.

Rigby, R. and Stasinopoulos, D. M (2013) Automatic smoothing parameter selection in GAMLSS with an application to centile estimation, \emph{Statistical methods in medical research}.

Stasinopoulos D. M. Rigby R.A. (2007) Generalized additive models for location scale and shape (GAMLSS) in R.
\emph{Journal of Statistical Software}, Vol. \bold{23}, Issue 7, Dec 2007, \url{http://www.jstatsoft.org/v23/i07}.

Stasinopoulos D. M., Rigby R.A., Heller G., Voudouris V., and De Bastiani F., (2017) \emph{Flexible Regression and Smoothing: Using GAMLSS in R}, Chapman and Hall/CRC. 

(see also \url{http://www.gamlss.org/}).

}
\author{Mikis Stasinopoulos, Vlasios Voudouris, Daniil Kiose}


\seealso{
\code{\link{gamlss}}
}
\examples{
m1 <- gamlss(R~pb(Fl)+pb(A), data=rent, family=GA)
# getting the Partial Efect function
pef <- getPEF(obj=m1,term="A", plot=TRUE)
# the value at 1980
pef(1980)
# the first derivative at 1980
pef(1980, deriv=1)
# the second derivative at 1980
pef(1980, deriv=2)
# plotting the first derivative
curve(pef(x, deriv=1), 1900,2000)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{regression}% use one of  RShowDoc("KEYWORDS")
