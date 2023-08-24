\name{edf}
\alias{edf}
\alias{edfAll}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Effective degrees of freedom from  gamlss model
}
\description{
The functions \code{edf()} and \code{edfAll()} can be used to obtained the effective degrees of freedom for 
different additive terms for the distribution parameters in a gamlss model. 
}
\usage{
edf(obj, what = c("mu", "sigma", "nu", "tau"),
    parameter= NULL, print = TRUE, ...)
edfAll(obj, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obj}{A gamlss fitted model}
  \item{what}{which of the four parameters \code{mu}, \code{sigma}, \code{nu} or \code{tau}.}
  \item{parameter}{equivalent to \code{what}}
  \item{print}{whether to print the label}
  \item{\dots}{for extra arguments}
}
\value{
The function \code{edfAll()} re turns a list of edf for  all the fitted parameters.
The function \code{edf()} a vector of edf.    
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
\author{
Mikis Stasinopoulos}
\note{The edf given are the ones fitted in the backfitting so the usually contained (depending on the additive term)
the contatnt and the linear part.
}



\seealso{
\code{\link{gamlss}}
}
\examples{
library(gamlss.data)
data(usair)
m1<- gamlss(y~pb(x1)+pb(x2)+pb(x6), data=usair)
edfAll(m1)
edf(m1)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{regression}
