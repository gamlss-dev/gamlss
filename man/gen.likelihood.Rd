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

Stasinopoulos D. M., Rigby R.A. and Akantziliotou C. (2006) Instructions on how to use the GAMLSS package in R.
Accompanying documentation in the current GAMLSS  help files, (see also  \url{http://www.gamlss.org/}). 

Stasinopoulos D. M. Rigby R.A. (2007) Generalized additive models for location scale and shape (GAMLSS) in R.
\emph{Journal of Statistical Software}, Vol. \bold{23}, Issue 7, Dec 2007, \url{http://www.jstatsoft.org/v23/i07}.

}
\author{
Mikis Stasinopoulos \email{d.stasinopoulos@londonmet.ac.uk}  Bob Rigby \email{r.rigby@londonmet.ac.uk} 
and Vlasios Voudouris}

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
