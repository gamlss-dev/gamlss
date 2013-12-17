\name{ra}
\alias{ra}

\title{Specify Simple Random Effect In A GAMLSS Formula}
\description{
This is an experimental smoother for use with factors in gamlss(). 
It allows the fitted values for a factor predictor to be shrunk towards the overall mean, 
where the amount of shrinking depends either on lambda, or on the equivalent degrees of freedom (df). 

This function is slightly more general, but considerably slower than the  \code{\link{random}} function .

}
\usage{
ra(xfactor, xvector = NULL, df = NULL, lambda = NULL, order = 0, 
   estimate = FALSE, expl = NULL, data1 = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{xfactor}{a factor defining the subjects grouping in a one factor random effect model term e.g. \code{xfactor=Subjects} }
  \item{xvector}{a variable if interaction with the \code{xfactor} is required \code{xvector} (experimental) }
  \item{df}{required equivalent degrees of freedom e.g. \code{df=10}}
  \item{lambda}{the smoothing parameter which is the reciprocal (i.e. inverse) of the variance of the random effect}
  \item{order}{the order of the difference in the matrix D, \code{order=1} is for simple random effects, 
               \code{order=2} is for random walk order 1 and \code{order=3} is for random walk order 2 }
  \item{estimate}{whether to estimate the lambda parameter within the backfitting iterations (very unreliable). Set by default to \code{estimate=FALSE}. 
                 [The lambda parameter can be more accurately estimated by selecting the corresponding smoothing degrees of freedom using \code{\link{find.hyper}}] }
  \item{expl}{this allows an explanatory variable at the subject level to be fitted e.g. \code{expl=~x1+x2}}
  \item{data1}{the data frame for the subject level variables \code{data1}}
}

\value{
  xfactor is returned with class "smooth", with an attribute named "call" which is to be evaluated in the backfitting  \code{additive.fit()} 
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
\author{ Mikis Stasinopoulos \email{d.stasinopoulos@londonmet.ac.uk}, Bob Rigby \email{r.rigby@londonmet.ac.uk}  }


\section{Warning }{This is experimental and likely to change soon}  
\seealso{ \code{\link{random}},  \code{\link{gamlss}}}
\examples{
data(aids)
attach(aids)
# fitting a loess curve with span=0.4 plus the a quarterly  effect 
aids1<-gamlss(y~lo(~x,span=0.4)+qrt,data=aids,family=PO) # 
# now we string the quarterly  effect using random 
aids2<-gamlss(y~lo(~x,span=0.4)+ra(qrt,df=2),data=aids,family=PO) # 
plot(x,y)
lines(x,fitted(aids1),col="red")
lines(x,fitted(aids2),col="purple")
rm(aids1,aids2)
detach(aids)
}
\keyword{regression}% 
