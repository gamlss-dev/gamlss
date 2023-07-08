\name{z.scores}
\alias{z.scores}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Z-scores for lms objects}
\description{
This creates z-scores for new values of y and x given a fitted \code{lms} object. 
}
\usage{
z.scores(object, y, x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{a \code{lms} fitted object}
  \item{y}{new y values}
  \item{x}{new x values}
}
\details{
This is simply a job that can be also done by \code{centiles.pred()}.
}
\value{
the required z-scores
}
\references{
Cole, T. J. (1994) Do growth chart centiles need a face lift?
 \emph{BMJ}, 308--641.

Cole,  T. J. and Green,  P. J. (1992) Smoothing reference centile curves: the LMS method and penalized likelihood, \emph{Statist. Med.} \bold{11}, 1305--1319

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
Mikis Stasinopoulos 
}

\seealso{
 \code{\link{centiles.pred}}
}
\examples{
\dontrun{
IND<-sample.int(7040, 1000, replace=FALSE)
db1 <- db[IND,]
plot(head~age, data=db1)
m0 <-  lms(head, age, data=db1,trans.x=TRUE )
z.scores(m0, x=c(2,15,30,40),y=c(45,50,56,63))}
}
\keyword{regression}

