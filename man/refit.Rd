\name{refit}
\alias{refit}

\title{Refit a GAMLSS model}
\description{
 This function refits a GAMLSS model. It is useful when the algorithm 
 has not converged after 20 outer iteration (the default value)}
\usage{
refit(object, ...)
}

\arguments{
  \item{object}{a GAMLSS fitted model which has not converged}
  \item{\dots}{for extra arguments}
}
\details{
This function is useful when the iterations have reach the maximum value set by the code(n.cyc) of the 
\code{gamlss.control} function and the model has not converged yet
}
\value{
 Returns a GAMLSS fitted model 
}
\references{ Rigby, R. A. and  Stasinopoulos D. M. (2005). Generalized additive models for location, scale and shape,(with discussion), 
\emph{Appl. Statist.}, \bold{54}, part 3, pp 507-554.

Stasinopoulos D. M., Rigby R.A. and Akantziliotou C. (2006) Instructions on how to use the GAMLSS package in R.
Accompanying documentation in the current GAMLSS  help files, (see also  \url{http://www.gamlss.org/}).

Stasinopoulos D. M. Rigby R.A. (2007) Generalized additive models for location scale and shape (GAMLSS) in R.
\emph{Journal of Statistical Software}, Vol. \bold{23}, Issue 7, Dec 2007, \url{http://www.jstatsoft.org/v23/i07}. 
 }
\author{ Mikis Stasinopoulos \email{d.stasinopoulos@londonmet.ac.uk}, Bob Rigby \email{r.rigby@londonmet.ac.uk} }
\note{The function \code{\link{update}} does a very similar job}


\seealso{ \code{\link{gamlss}}, \code{\link{update.gamlss}} }
\examples{
data(aids)
h<-gamlss(y~poly(x,3)+qrt, family=PO, data=aids) # 
refit(h)
rm(h)
}
\keyword{regression}% 
