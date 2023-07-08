\name{gamlss-package}
\alias{gamlss-package}
\docType{package}
\title{\packageTitle{gamlss}}
\description{
\packageDescription{gamlss}
}
\details{
The DESCRIPTION file:
\packageDESCRIPTION{gamlss}
\packageIndices{gamlss}
}
\author{
\packageAuthor{gamlss}

Maintainer: \packageMaintainer{gamlss}
}
\references{
Nelder, J. A. and Wedderburn, R. W. M. (1972). Generalized
linear models. \emph{J. R. Statist. Soc.} A., \bold{135} 370-384.

Hastie, T. J. and Tibshirani, R. J. (1990). \emph{Generalized Additive
Models}. Chapman and Hall, London.

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

\keyword{ package }
\keyword{ distribution }
\keyword{ regression }

\seealso{
\code{\link[gamlss.dist]{gamlss.dist}} 
}
\examples{
data(abdom)
mod<-gamlss(y~pb(x),sigma.fo=~pb(x),family=BCT, data=abdom, method=mixed(1,20))
plot(mod)
rm(mod)
}
