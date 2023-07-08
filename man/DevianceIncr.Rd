\name{devianceIncr}
\alias{devianceIncr}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{The global deviance increment}
\description{ The  global deviance increment is the contribution of each individual observation to the global deviance.   
The function \code{devianceIncr()} can be used to extract the global deviance increment for a fitted gamlss model  or for a new (test/validation) data set. Large values for global deviance increment indicate a bad fit and for new data a bad prediction. 
}
\usage{
devianceIncr(obj, newdata = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obj}{a gamlss object}
  \item{newdata}{test data set to check the global deviance increment.}
}

\value{
Returns a vector of the global deviance increments for each observation.
}
\references{
Rigby, R. A. and Stasinopoulos D. M. (2005). Generalized additive models for location, scale and shape,(with discussion), Appl. Statist., 54, part 3, 1-38.

Rigby, R. A., Stasinopoulos, D. M.,  Heller, G. Z.,  and De Bastiani, F. (2019)
	\emph{Distributions for modeling location, scale, and shape: Using GAMLSS in R}, Chapman and Hall/CRC. An older version can be found in \url{https://www.gamlss.com/}.

Stasinopoulos D. M. Rigby R.A. (2007) Generalized additive models for location scale and shape (GAMLSS) in R.
\emph{Journal of Statistical Software}, Vol. \bold{23}, Issue 7, Dec 2007, \url{https://www.jstatsoft.org/v23/i07/}.

Stasinopoulos D. M., Rigby R.A., Heller G., Voudouris V., and De Bastiani F., (2017)
\emph{Flexible Regression and Smoothing: Using GAMLSS in R},  Chapman and Hall/CRC.  

(see also \url{https://www.gamlss.com/}).
}
\author{Mikis Stasinopoulos}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{\link{deviance}}
}
\examples{
#-----------------------------------------------------------------
# Count data set
# fit Poisson model
h1 <- gamlss(Claims~L_Population+L_Accidents+L_KI+L_Popdensity, 
             data=LGAclaims, family=PO)
p1<-devianceIncr(h1)
# fit negative binomial model
h2 <- gamlss(Claims~L_Population+L_Accidents+L_KI+L_Popdensity, 
             data=LGAclaims, family=NBI)
p2<-devianceIncr(h2)
# comparing using boxplots
boxplot(cbind(p1,p2))
# comparing using emphirical cdf
plot(ecdf(p1))
lines(ecdf(p2), col=2)
# comparing agaist the y-values
plot(p1~LGAclaims$Claims, pch=20, col="gray")
points(p2~LGAclaims$Claims, pch="-", col="orange")
#----------------------------------------------------------------
# Continuous data sets
\dontrun{
m1 <- gamlss(head~pb(age), data=db[1:6000,])
p1<-devianceIncr(m1)
m2 <- gamlss(head~pb(age), sigma.fo=~pb(age), nu.fo=~pb(age), 
      tau.fo=~pb(age), data=db[1:6000,], family=BCT)
p2<-devianceIncr(m2)
# comparing using summaries
summary(p1); summary(p2)
# comparing using boxplots
boxplot(cbind(p1,p2))
# comparing using histograms
hist(p1, col=rgb(1,0,0,0.5), xlim=c(0,50), breaks=seq(0,50,2))
hist(p2, col=rgb(0,0,1,0.5), add=T)
# comparing using emphirical cdf
plot(ecdf(p1))
lines(ecdf(p2), col=2)
}
#----------------------------------------------------------------
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{regression}% use one of  RShowDoc("KEYWORDS")

