\name{gen.likelihood}
\alias{gen.likelihood}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
A function to generate the likelihood function from a GAMLSS object
}
\description{
This function generate a function with argument the parameters of the GAMLSS model which can evaluate the log-likelihood function. 
}
\usage{
gen.likelihood(object)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{A gamlss fitted model}
}
\details{
The purpose of this function is to help the function vcov() to get he right Hessian matrix after a model has fitted. Note that at the momment smoothing terms are consideted as fixed. 
}
\value{
A function of the log-likelihood
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
Mikis Stasinopoulos, Bob Rigby and Vlasios Voudouris}

\seealso{
\code{\link{vcov}}
}
\examples{
data(aids)
m1 <- gamlss(y~x+qrt, data=aids, family=NBI)
logL<-gen.likelihood(m1)
logL()
logLik(m1)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{regression}
