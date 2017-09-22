\name{LR.test}
\alias{LR.test}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Likelihood Ratio test for nested GAMLSS models
}
\description{
The function performs  a likelihood ration test for two nested fitted model.  
}
\usage{
LR.test(null, alternative, print = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{null}{The null hypothesis (simpler) fitted model 
}
  \item{alternative}{
The alternative hypothesis (more complex) fitted model
}
  \item{print}{whether to print or save the result}
} \details{ Warning: no checking whether the models are nested 
is performed. } \value{ If \code{print=FALSE} a list with 
\code{chi}, \code{df} and \code{p.val} is produced.
}
\references{
Rigby, R. A. and  Stasinopoulos D. M. (2005). Generalized additive models for location, scale and shape,(with discussion), 
\emph{Appl. Statist.}, \bold{54}, part 3, pp 507-554.

Stasinopoulos D. M. Rigby R.A. (2007) Generalized additive models for location scale and shape (GAMLSS) in R.
\emph{Journal of Statistical Software}, Vol. \bold{23}, Issue 7, Dec 2007, \url{http://www.jstatsoft.org/v23/i07}.

Stasinopoulos D. M., Rigby R.A., Heller G., Voudouris V., and De Bastiani F., (2017) \emph{Flexible Regression and Smoothing: Using GAMLSS in R}, Chapman and Hall/CRC. 

(see also \url{http://www.gamlss.org/}).
}
\author{
Mikis Stasinopoulos \email{mikis.stasinopoulos@gamlss.org}}


\seealso{
\code{\link{gamlss}}, \code{\link{dropterm} }
}
\examples{
data(usair)
m0<-gamlss(y~x1+x2, data=usair)
m1<-gamlss(y~x1+x2+x3+x4, data=usair)
LR.test(m0,m1)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{regression}
