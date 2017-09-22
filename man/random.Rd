\name{random}
\alias{random}
\alias{re}

\title{Specify a random intercept model in a GAMLSS formula}

\description{They are two functions for fitting random effects wthin a  GAMLSS model, \code{random()} and \code{re()}. 

The function \code{random()} is based on the original \code{random()} function of Trevor Hastie in the package \code{gam}. TIn our version the function has been modified to allow a "local" maximum likelihood estimation of the smoothing parameter \code{lambda}. This method is  equivalent to the PQL method of Breslow and Clayton (1993) applied at the local iterations of the algorithm. In fact for a GLM model and a simple random effect it is equivalent to \code{glmmPQL()} function in the package \code{MASS} see Venables and Ripley (2002).  Venables and Ripley (2002) claimed that this iterative method was first introduced by Schall (1991). Note that in order for the "local" maximum likelhood estimation procedure to operate both argument \code{df} and \code{lambda} has to be \code{NULL}.

The function \code{re()} is an interface for calling the \code{lme()} function of the package  \pkg{nlme}. This gives the user the abilty to fit comlpicated random effect models while the assumtion of the normal distribution for the response variable is relaxed.  The theoretical justification cames again from the fact that this is a PQL method,  Breslow and Clayton (1993). 

}
\usage{
random(x, df = NULL, lambda = NULL, start=10)

re(fixed = ~1, random = NULL, correlation = NULL, method = "ML", ...)
}

\arguments{
  \item{x}{a factor }
  \item{df}{the target degrees of freedom}
  \item{lambda}{the smoothing parameter lambda which can be viewed as a shrinkage parameter.}
  \item{start}{starting value for lambda if local Maximul likelihood is used.}
  \item{fixed}{a formula specify the fixed effects of the \code{lme()} model. This, in most cases  can be also included in the \code{gamlss} parameter formula}
  \item{random}{a formula or list specufying the random effect  part of the model as in  \code{lme()} function}
  \item{correlation}{the correlation structure of the \code{lme()} model}
  \item{method}{which method, "ML" (the default), or "REML"}
  \item{\dots}{this can be used to pass arguments for \code{lmeControl()}
}
  }

\details{
The function \code{random()} can be seen as  a smoother for use with factors in gamlss(). 
It allows the fitted values for a factor predictor to be shrunk towards the overall mean, 
where the amount of shrinking depends either on lambda, or on the equivalent degrees of freedom or on the estimated sigma parameter (default).  Similar in spirit to smoothing splines, this fitting method can be justified on Bayesian grounds or by a random effects model. Note that the behavier of the function is different from the original Hastie function. Here the function behaves as follows: i) if both \code{df} and \code{lambda} are \code{NULL} then the PQL method is used
ii) if  \code{lambda} is not \code{NULL},  \code{lambda} is used for fitting
iii) if  \code{lambda} is  \code{NULL} and   \code{df} is not \code{NULL} then \code{df}
 is used for fitting. 
 
Since factors are coded by model.matrix() into a set of contrasts, care has been taken to add an appropriate "contrast" 
attribute to the output of random(). This zero contrast results in a column of zeros in the model matrix,  which is aliased with any column and is hence ignored.


The use of the function \code{re()} requires knowledge of the use of the function \code{lme()} of the package \pkg{nlme} for the specification of the appropriate random effect model. Some care should betaken whether the data set is  
}
\value{
 x is returned with class "smooth", with an attribute named "call" which is to be evaluated in the backfitting  \code{additive.fit()} 
   called by \code{gamlss()}
}
\references{
Breslow, N. E. and Clayton, D. G. (1993) Approximate inference in generalized linear mixed models. \emph{Journal of the American Statistical Association} \bold{88}, 9???25.

Chambers, J. M. and Hastie, T. J. (1991). \emph{Statistical Models in S}, Chapman and Hall, London. 

Pinheiro, Jose C and Bates, Douglas M (2000) \emph{Mixed effects models in S and S-PLUS}
Springer.

Rigby, R. A. and  Stasinopoulos D. M. (2005). Generalized additive models for location, scale and shape,(with discussion), 
\emph{Appl. Statist.}, \bold{54}, part 3, pp 507-554.

Schall, R. (1991) Estimation in generalized linear models with random effects. \emph{Biometrika} \bold{78}, 719???727.

Stasinopoulos D. M. Rigby R.A. (2007) Generalized additive models for location scale and shape (GAMLSS) in R.
\emph{Journal of Statistical Software}, Vol. \bold{23}, Issue 7, Dec 2007, \url{http://www.jstatsoft.org/v23/i07}.

Stasinopoulos D. M., Rigby R.A., Heller G., Voudouris V., and De Bastiani F., (2017) \emph{Flexible Regression and Smoothing: Using GAMLSS in R}, Chapman and Hall/CRC. (see also \url{http://www.gamlss.org/}).

Venables, W. N. and Ripley, B. D. (2002) \emph{Modern Applied Statistics with S}. Fourth edition. Springer.
}
            
\author{For \code{re()} Mikis Stasinopoulos and Marco Enea and for \code{random()} Trevor Hastie (amended by Mikis Stasinopoulos), }

\seealso{\code{\link{gamlss}}, \code{\link{gamlss.random}}}

\examples{
#------------- Example 1 from Pinheiro and Bates (2000) page 15-----------------
# bring nlme
library(nlme)
data(ergoStool)
# lme model
l1<-lme(effort~Type, data=ergoStool, random=~1|Subject, method="ML")
# use random() 
t1<-gamlss(effort~Type+random(Subject), data=ergoStool )
# use re() with fixed effect within re()
t2<-gamlss(effort~re(fixed=~Type, random=~1|Subject), data=ergoStool )
# use re() with fixed effect in gamlss formula
t3<-gamlss(effort~Type+re(random=~1|Subject), data=ergoStool )
# compare lme fitted values with random
plot(fitted(l1), fitted(t1))
# compare lme fitted values with random
plot(fitted(l1), fitted(t2))
lines(fitted(l1), fitted(t3), col=2)
# getting the fitted coefficients 
getSmo(t2)
#-------------------------------------------------------------------------------
\dontrun{
#-------------Example 2 Hodges data---------------------------------------------
data(hodges)
plot(prind~state, data=hodges)
m1<- gamlss(prind~random(state), sigma.fo=~random(state), nu.fo=~random(state), 
            tau.fo=~random(state), family=BCT, data=hodges)
m2<- gamlss(prind~re(random=~1|state), sigma.fo=~re(random=~1|state), 
            nu.fo=~re(random=~1|state), tau.fo=~re(random=~1|state), family=BCT, 
            data=hodges)
# comparing the fitted effective degrees of freedom
m1$mu.df
m2$mu.df
m1$sigma.df
m2$sigma.df
m1$nu.df
m2$nu.df
m1$tau.df
m2$tau.df
# random effect for tau is not needed
m3<- gamlss(prind~random(state), sigma.fo=~random(state), nu.fo=~random(state),  
            family=BCT, data=hodges, start.from=m1)
plot(m3)
# term plots work for random but not at the moment for re()
op <- par(mfrow=c(2,2))
term.plot(m3, se=TRUE)
term.plot(m3, se=TRUE, what="sigma")
term.plot(m3, se=TRUE, what="nu")
par(op)
# getting information from a fitted lme object
coef(getSmo(m2))
ranef(getSmo(m2))
VarCorr(getSmo(m2))
summary(getSmo(m2))
intervals(getSmo(m2))
fitted(getSmo(m2))
fixef(getSmo(m2))
#  plotting 
plot(getSmo(m2))
qqnorm(getSmo(m2))
#----------------Example 3 from Pinheiro and Bates (2000) page 42---------------
data(Pixel)
l1 <- lme(pixel~ day+I(day^2), data=Pixel, random=list(Dog=~day, Side=~1),
          method="ML")
# this will fail 
#t1<-gamlss(pixel~re(fixed=~day+I(day^2), random=list(Dog=~day, Side=~1)), 
#           data=Pixel)
# but this  is working 
t1<-gamlss(pixel~re(fixed=~day+I(day^2), random=list(Dog=~day, Side=~1), 
                    opt="optim"), data=Pixel)
plot(fitted(l1)~fitted(t1))
#---------------Example 4 from Pinheiro and Bates (2000)page 146----------------
data(Orthodont)
l1 <- lme(distance~ I(age-11), data=Orthodont, random=~I(age-11)|Subject,
           method="ML")

t1<-gamlss(distance~I(age-11)+re(random=~I(age-11)|Subject), data=Orthodont)
plot(fitted(l1)~fitted(t1))
# checking the model
plot(t1)
wp(t1, ylim.all=2)
# two observation fat try LO
t2<-gamlss(distance~I(age-11)+re(random=~I(age-11)|Subject,  opt="optim", 
     numIter=100), data=Orthodont, family=LO)
plot(t2)
wp(t2,ylim.all=2)
# a bit better but not satisfactory Note that  3 paramters distibutions fail
#------------example 5 from Venable and Ripley (2002)--------------------------
library(MASS)
data(bacteria)
summary(glmmPQL(y ~ trt + I(week > 2), random = ~ 1 | ID,
                family = binomial, data = bacteria))
s1 <- gamlss(y ~ trt + I(week > 2)+random(ID), family = BI, data = bacteria)
s2 <- gamlss(y ~ trt + I(week > 2)+re(random=~1|ID), family = BI, 
             data = bacteria)
s3 <- gamlss(y ~ trt + I(week > 2)+re(random=~1|ID, method="REML"), family = BI, 
             data = bacteria)
# the esimate of the random effect sd sigma_b 
sqrt(getSmo(s1)$tau2)
getSmo(s2)
getSmo(s3)
#-------------Example 6 from Pinheiro and Bates (2000) page 239-244-------------
# using corAR1()
data(Ovary)
# AR1 
l1 <- lme(follicles~sin(2*pi*Time)+cos(2*pi*Time), data=Ovary, 
          random=pdDiag(~sin(2*pi*Time)), correlation=corAR1())
# ARMA
l2 <- lme(follicles~sin(2*pi*Time)+cos(2*pi*Time), data=Ovary, 
          random=pdDiag(~sin(2*pi*Time)), correlation=corARMA(q=2))
# now gamlss
# AR1 
t1 <- gamlss(follicles~re(fixed=~sin(2*pi*Time)+cos(2*pi*Time), 
                         random=pdDiag(~sin(2*pi*Time)),
                         correlation=corAR1()), data=Ovary)
plot(fitted(l1)~fitted(t1))
# ARMA
t2 <- gamlss(follicles~re(fixed=~sin(2*pi*Time)+cos(2*pi*Time), 
                          random=pdDiag(~sin(2*pi*Time)),
                          correlation=corARMA(q=2)), data=Ovary)
plot(fitted(l2)~fitted(t2))
AIC(t1,t2)
wp(t2, ylim.all=1)
#-------------------------------------------------------------------------------  
}
}
\keyword{regression}% 
