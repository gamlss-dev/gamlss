\name{print.gamlss}
\alias{print.gamlss}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Prints a GAMLSS fitted model }
\description{
\code{print.gamlss} is the GAMLSS specific method for the generic function \code{print} which prints 
      objects returned by modelling functions. 
}
\usage{
\method{print}{gamlss}(x, digits = max(3, getOption("digits") - 3), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{a GAMLSS fitted model}
  \item{digits}{the number of significant digits to use when printing }
  \item{\dots}{for extra arguments}
}

\value{Prints a gamlss object}
\references{ Rigby, R. A. and  Stasinopoulos D. M. (2005). Generalized additive models for location, scale and shape,(with discussion), 
\emph{Appl. Statist.}, \bold{54}, part 3, pp 507-554.

Stasinopoulos D. M., Rigby R.A. and Akantziliotou C. (2006) Instructions on how to use the GAMLSS package in R.
Accompanying documentation in the current GAMLSS  help files,(see also  \url{http://www.gamlss.org/}). 

Stasinopoulos D. M. Rigby R.A. (2007) Generalized additive models for location scale and shape (GAMLSS) in R.
\emph{Journal of Statistical Software}, Vol. \bold{23}, Issue 7, Dec 2007, \url{http://www.jstatsoft.org/v23/i07}.
}
\author{Mikis Stasinopoulos \email{d.stasinopoulos@londonmet.ac.uk}, Bob Rigby \email{r.rigby@londonmet.ac.uk} and Calliope Akantziliotou }

\seealso{  \code{\link{gamlss}}, \code{\link{deviance.gamlss}},  \code{\link{fitted.gamlss}} }
\examples{
data(aids)
h<-gamlss(y~poly(x,3)+qrt, family=PO, data=aids) 
print(h) # or just h
rm(h)
}
\keyword{regression}% 
