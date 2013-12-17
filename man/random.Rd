\name{random}
\alias{random}

\title{Specify a simple random effect in a GAMLSS Formula}
\description{ Includes random effect terms in an GAMLSS model. The function is based on the original  \code{random()} function of Trevor Hastie in the package \code{gam}. This version of the function has been modified to allow a "local" maximum likelihood estimation of the smoothing parameter \code{lambda}. This method is  equivalent to the PQL method of Breslow and Clayton (1993) applied at the local iterations of the algorithm. In fact for a GLM model and a simple random effect it is equivalent to \code{glmmPQL()} function in the package \code{MASS} see Venables and Ripley (2002).  Venables and Ripley (2002) claimed that this iterative method was first introduced by Schall (1991). Note that in order for the "local" maximum likelhood estimation procedure to operate both argument \code{df} and \code{lambda} has to be \code{NULL}.}

\usage{random(x, df = NULL, lambda = NULL, start=10)}

\arguments{
  \item{x}{a factor }
  \item{df}{the target degrees of freedom}
  \item{lambda}{the smoothing parameter lambda which can be viewed as a shrinkage parameter.}
  \item{start}{starting value for lambda if local Maximul likelihood is used.}
  }

\details{
This is a smoother for use with factors in gamlss(). 
It allows the fitted values for a factor predictor to be shrunk towards the overall mean, 
where the amount of shrinking depends either on lambda, or on the equivalent degrees of freedom. 
Similar in spirit to smoothing splines, this fitting method can be justified on Bayesian grounds or by a random effects model.

Note that the behavier of the function is different from the original Hastie function. Here the function behaves as follows: i) if both \code{df} and \code{lambda} are \code{NULL} then the PQL method is used
ii) if  \code{lambda} is not \code{NULL},  \code{lambda} is used for fitting
iii) if  \code{lambda} is  \code{NULL} and   \code{df} is not \code{NULL} then \code{df}
 is used for fitting. 
 
Since factors are coded by model.matrix() into a set of contrasts, care has been taken to add an appropriate "contrast" 
attribute to the output of random(). This zero contrast results in a column of zeros in the model matrix, 
which is aliased with any column and is hence ignored
}
\value{
 x is returned with class "smooth", with an attribute named "call" which is to be evaluated in the backfitting  \code{additive.fit()} 
   called by \code{gamlss()}
}
\references{
Breslow, N. E. and Clayton, D. G. (1993) Approximate inference in generalized linear mixed models. \emph{Journal of the American Statistical Association} \bold{88}, 9???25.

Chambers, J. M. and Hastie, T. J. (1991). \emph{Statistical Models in S}, Chapman and Hall, London. 

Rigby, R. A. and  Stasinopoulos D. M. (2005). Generalized additive models for location, scale and shape,(with discussion), 
\emph{Appl. Statist.}, \bold{54}, part 3, pp 507-554.

Schall, R. (1991) Estimation in generalized linear models with random effects. \emph{Biometrika} \bold{78}, 719???727.

Stasinopoulos D. M., Rigby R.A. and Akantziliotou C. (2006) Instructions on how to use the GAMLSS package in R.
Accompanying documentation in the current GAMLSS  help files, (see also  \url{http://www.gamlss.org/}).

Stasinopoulos D. M. Rigby R.A. (2007) Generalized additive models for location scale and shape (GAMLSS) in R.
\emph{Journal of Statistical Software}, Vol. \bold{23}, Issue 7, Dec 2007, \url{http://www.jstatsoft.org/v23/i07}.

Venables, W. N. and Ripley, B. D. (2002) \emph{Modern Applied Statistics with S}. Fourth edition. Springer.
}
            
\author{ Trevor Hastie (amended by Mikis Stasinopoulos)}


\seealso{\code{\link{gamlss}}, \code{\link{gamlss.random}}}

\examples{
\dontrun{
data(hodges)
plot(prind~state, data=hodges)
m1<- gamlss(prind~random(state), sigma.fo=~random(state), nu.fo=~random(state), 
             tau.fo=~random(state), family=BCT, data=hodges)
edfAll(m1)
# radnom effect for tau is not needed
m2<- gamlss(prind~random(state), sigma.fo=~random(state), nu.fo=~random(state),  
            family=BCT, data=hodges, start.from=m1)
edfAll(m2)
plot(m2)
#op<-par(mfrow=c(3,1))
#term.plot(m2, se=TRUE)
#term.plot(m2, se=TRUE, what="sigma")
#term.plot(m2, se=TRUE, what="nu")
#par(op)
# the example from Venable and Ripley (2002)
library(MASS)
data(bacteria)
library(nlme)
summary(glmmPQL(y ~ trt + I(week > 2), random = ~ 1 | ID,
                family = binomial, data = bacteria))
s1 <- gamlss(y ~ trt + I(week > 2)+random(ID), family = BI, data = bacteria)
# the esimate of sigma 
sqrt(s1$mu.coefSmo[[1]]$sig2)
# the esimate of random effect  sigma
sqrt(s1$mu.coefSmo[[1]]$tau2)
}
}
\keyword{regression}% 
