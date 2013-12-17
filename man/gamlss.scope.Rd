\name{gamlss.scope}
\alias{gamlss.scope}


\title{Generate a Scope Argument for Stepwise GAMLSS}
\description{
  Generates a scope argument for a stepwise GAMLSS.
}
\usage{
gamlss.scope(frame, response = 1, smoother = "cs", arg = NULL, form = TRUE)
}

\arguments{
  \item{frame}{a data or model frame}
  \item{response}{ which variable is the response; the default is the first}
  \item{smoother}{what smoother to use; default is \code{cs}}
  \item{arg}{any additional arguments required by the smoother}
  \item{form}{should a formula be returned (default), or else a character version of the formula }
}
\details{
 Each formula describes an ordered regimen of terms, each of which is eligible on their own for 
 inclusion in the gam model. 
 One of the terms is selected from each formula by step.gam. 
 If a 1 is selected, that term is omitted. 
}
\value{
a list of formulas is returned, one for each column in frame (excluding the response). 
For a numeric variable, say x1, the formula is

~ 1 + x1 + cs(x1)

If x1 is a factor, the last smooth term is omitted.
}
\references{Chambers, J. M. and Hastie, T. J. (1991). \emph{Statistical Models in S}, Chapman and Hall, London. 

Rigby, R. A. and  Stasinopoulos D. M. (2005). Generalized additive models for location, scale and shape,(with discussion), 
\emph{Appl. Statist.}, \bold{54}, part 3, pp 507-554.

Stasinopoulos D. M., Rigby R.A. and Akantziliotou C. (2006) Instructions on how to use the GAMLSS package in R.
Accompanying documentation in the current GAMLSS  help files, (see also  \url{http://www.gamlss.org/}).

Stasinopoulos D. M. Rigby R.A. (2007) Generalized additive models for location scale and shape (GAMLSS) in R.
\emph{Journal of Statistical Software}, Vol. \bold{23}, Issue 7, Dec 2007, \url{http://www.jstatsoft.org/v23/i07}.
}
\author{Mikis Stasinopoulos: a modified function from Statistical Models in S}


\seealso{  \code{\link{stepGAIC}}}
\examples{
data(usair)
gs1<-gamlss.scope(model.frame(y~x1+x2+x3+x4+x5+x6, data=usair))
gs2<-gamlss.scope(model.frame(usair))
gs1
gs2
gs3<-gamlss.scope(model.frame(usair), smooth="fp", arg="3")
gs3 
}
\keyword{regression}
