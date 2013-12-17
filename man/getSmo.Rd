\name{getSmo}
\alias{getSmo}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Extracting Smoother information from a GAMLSS fitted object
}
\description{
The function \code{getSmo()} extracts information from a fitted smoothing additive term.  
}
\usage{
getSmo(object, what = c("mu", "sigma", "nu", "tau"), 
       which = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
a GAMLSS fitted model
}
  \item{what}{
which distribution parameter is required, default what="mu"
}
  \item{which}{
which smoothing term
}
}
\details{
This function facilitates the extraction of information from a fitted additive terms. For example \code{getSmo(m1,"sigma",2)} is equivalent of \code{m1$sigma.coefSmo[[2]]}. To get the actual fitted values type \code{m1$sigma.s[[2]]}
}
\value{
A list containing information about a fitted smoother or a fitted objects 
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
Mikis Stasinopoulos and Bob Rigby
}


%% ~Make other sections like Warning with \section{Warning }{....} ~


\examples{
data(usair)
t1<-gamlss(y~x1+pb(x5)+pb(x6), data=usair, family=GA)
# get the value for lambda for the second fitted term in mu
getSmo(t1, "mu", 2)$lambda
  }
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{regression}

