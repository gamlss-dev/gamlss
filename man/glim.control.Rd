\name{glim.control}
\alias{glim.control}

\title{ Auxiliary for Controlling the inner algorithm in a GAMLSS Fitting }
\description{
  Auxiliary function used for the inner iteration of \code{gamlss} algorithm. Typically
     only used when calling \code{gamlss} function through the option \code{i.control}. 
}
\usage{
glim.control(cc = 0.001, cyc = 50,  glm.trace = FALSE, 
             bf.cyc = 30, bf.tol = 0.001, bf.trace = FALSE, 
             ...)
}

\arguments{
  \item{cc}{the convergence criterion for the algorithm}
  \item{cyc}{the number of cycles of the algorithm}
  \item{glm.trace}{whether to print at each iteration (TRUE) or not (FALSE) }
  \item{bf.cyc}{the number of cycles of the backfitting algorithm}
  \item{bf.tol}{the convergence criterion (tolerance level) for the backfitting algorithm }
  \item{bf.trace}{whether to print at each iteration (TRUE) or not (FALSE, the default)}
  \item{\dots}{for extra arguments}
}


\value{A list with the arguments as components}

\references{ 
Rigby, R. A. and  Stasinopoulos D. M. (2005). Generalized additive models for location, scale and shape, (with discussion), 
\emph{Appl. Statist.}, \bold{54}, part 3, pp 507-554.

Stasinopoulos D. M., Rigby R.A. and Akantziliotou C. (2006) Instructions on how to use the GAMLSS package in R.
Accompanying documentation in the current GAMLSS  help files, (see also  \url{http://www.gamlss.org/}).  

Stasinopoulos D. M. Rigby R.A. (2007) Generalized additive models for location scale and shape (GAMLSS) in R.
\emph{Journal of Statistical Software}, Vol. \bold{23}, Issue 7, Dec 2007, \url{http://www.jstatsoft.org/v23/i07}.
}

\author{Mikis Stasinopoulos \email{d.stasinopoulos@londonmet.ac.uk}, Bob Rigby \email{r.rigby@londonmet.ac.uk}}
 
\seealso{ \code{\link{gamlss}} }
\examples{
data(aids)
con<-glim.control(glm.trace=TRUE)
h<-gamlss(y~poly(x,3)+qrt, family=PO, data=aids, i.control=con) # 
rm(h,con)

}
\keyword{regression}% 
