\name{ridge}
\alias{ridge}
\alias{ri}

\title{Specify Ridge regression  In A GAMLSS Formula}
\description{
The function \code{ri()} allow the user to fit a ridge regression within GAMLSS. 
It allows the coefficients of a set of explanatory variables to be shrunk towards an overall zero, 
where the amount of shrinking depends either on lambda, or on the equivalent degrees of freedom (df). 
The function \code{ri()} can estimate lambda a local REML-algorithm.
}
\usage{
ridge(X, df = NULL, lambda = NULL, order = 0)
ri(X, df = NULL, lambda = NULL, order = 0, start = 10)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{ A standardized (mean=0, sd-1)  matrix  \code{X} of explanatory variables  }
  \item{df}{the effective degrees of freedom \code{df}  }
  \item{lambda}{the smoothing parameter \code{lambda}}
  \item{order}{ the \code{order} of the difference applied to the coefficients }
   \item{start}{the lambda starting value if the Schall's EM-algorithm is used}
}

\value{
  x is returned with class "smooth", with an attribute named "call" which is to be evaluated in the backfitting  \code{additive.fit()} 
   called by \code{gamlss()}
}
\references{
Rigby, R. A. and  Stasinopoulos D. M. (2005). Generalized additive models for location, scale and shape,(with discussion), 
\emph{Appl. Statist.}, \bold{54}, part 3, pp 507-554.

Stasinopoulos D. M., Rigby R.A. and Akantziliotou C. (2006) Instructions on how to use the GAMLSS package in R.
Accompanying documentation in the current GAMLSS  help files, (see also  \url{http://www.gamlss.org/}).
 
Stasinopoulos D. M. Rigby R.A. (2007) Generalized additive models for location scale and shape (GAMLSS) in R.
\emph{Journal of Statistical Software}, Vol. \bold{23}, Issue 7, Dec 2007, \url{http://www.jstatsoft.org/v23/i07}. 
 }
\author{ Mikis Stasinopoulos \email{d.stasinopoulos@londonmet.ac.uk}, Bob Rigby \email{r.rigby@londonmet.ac.uk} and Paul Eilers}


\seealso{ \code{\link{ra}},  \code{\link{gamlss}}}
\examples{
data(usair)
# standarized the x's
for (i in 2:7) usair[,i] <- (usair[,i]-mean(usair[,i]))/sd(usair[,i])
# create a matrix 
usair$X <-as.matrix(usair[,c("x1","x2", "x3","x4", "x5", "x6")])
m1<- gamlss(y~ridge(X, df=3), data=usair, family=GA)
m2<- gamlss(y~ridge(X, lambda=10), data=usair, family=GA)
# plotting the coeficients as function of the df 
df <- seq(0,6,0.5)
resmat<-matrix(0, nrow=length(df), ncol=6)
for (i in 1:length(df)) 
{
resmat[i,] <- gamlss(y~ridge(X, df=df[i]), data=usair)$mu.coefSmo[[1]][["coef"]]
}
colnames(resmat)<-colnames(usair$X)
plot(1:length(df), type="n", xlim=c(-.5,6.5), ylim=c(-23, 38))
for (i in 1:6)
{
lines(resmat[,i]~df, col="purple")
points(resmat[,i]~df)
}
lines(rep(0,length(df))~df, col="black")
# now estimating the lambda
 g1<-gamlss(y~ri(X), data=usair)
  m1$mu.coefSmo
}
\keyword{regression}% 
