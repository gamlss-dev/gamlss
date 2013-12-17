\name{stepGAIC}

\alias{stepGAIC}
\alias{stepGAIC.CH}
\alias{stepGAIC.VR}
\alias{stepGAICAll.A}
\alias{stepGAICAll.B}

\alias{stepTGD}

\title{Choose a model by GAIC in a Stepwise Algorithm}

\description{
The function \code{stepGAIC()} performs stepwise 
model selection using a Generalized Akaike Information Criterion.
The function \code{stepGAIC()} calls one of the two functions \code{stepGAIC.VR()} or \code{stepGAIC.CH()} depending on the argument 
\code{additive}. 
The function \code{stepGAIC.VR()} is based on the function \code{stepAIC()} given in the library MASS of Venables and Ripley (2002).
The function  \code{stepGAIC.CH} is based on the S function \code{step.gam()} (see Chambers and Hastie (1991)) 
and it is more suited for model with smoothing additive terms, (see below comments for the additive function \code{pb()}).  
Both functions have been adapted to work with  gamlss objects.  The main difference for the user is the \code{scope} 
argument, see below. 
If the \code{stepGAIC()} is called with the argument \code{additive=FALSE} then the  \code{stepGAIC.VR()} is called else the 
\code{stepGAIC.CH()}.

While the functions \code{stepGAIC.VR()} and \code{stepGAIC.CH()} are used to build models for individual parameters of the distribution 
of the response variable, the functions \code{stepGAICAll.A()} and \code{stepGAICAll.A()}  are building a model for all 
the parameters.
Both the functions  \code{stepGAICAll.A()} and \code{stepGAICAll.B()} are based on \code{stepGAIC.VR()}.
The use two different  strategies for selecting a appropriate final model.    \code{stepGAICAll.A()} has   the following
strategy: 
 
Stategy A:

 
i) build a model for \code{mu} using a forward approach. 

ii) given the model for \code{mu} build a model for \code{sigma}  (forward)

iii) given the models for  \code{mu} and \code{sigma} build a model for \code{nu} (forward) 

iv)  given the models for  \code{mu}, \code{sigma} and \code{nu} build a model for \code{tau} (forward) 

v) given the models for  \code{mu}, \code{sigma},  \code{nu} and \code{tau} check whether the terms for \code{nu} 
are needed using backward elimination. 

vi) given the models for  \code{mu}, \code{sigma},  \code{nu} and \code{tau} check whether the terms for \code{sigma} 
are needed (backward).

vii) given the models for  \code{mu}, \code{sigma},  \code{nu} and \code{tau} check whether the terms for \code{mu} 
are needed (backward).

Note for this strategy to work the \code{scope} argument should be set appropriately. 


\code{stepGAICAll.B()} uses the same procedure as the function  \code{stepGAIC.VR()} but each term in the scope is fitted 
to ALL the parameters of the distribution, rather than the one specified  by the argument \code{what} of \code{stepGAIC.VR()}.
 

}
\usage{

stepGAIC.VR(object, scope, direction = c("both", "backward", "forward"), 
         trace = T, keep = NULL, steps = 1000, scale = 0, 
         what = c("mu", "sigma", "nu", "tau"), k = 2, ...)

stepGAIC.CH(object, scope = gamlss.scope(model.frame(object)), 
            direction = c("both", "backward", "forward"), trace = T, keep = NULL, 
            steps = 1000, what = c("mu", "sigma", "nu", "tau"), k = 2, ...)


stepGAIC(object, scope = gamlss.scope(model.frame(object)), 
          direction = c("both", "backward", "forward"), 
          trace = T, keep = NULL, steps = 1000, 
          what = c("mu", "sigma", "nu", "tau"), k = 2, 
          additive = FALSE, ...)

stepGAICAll.A(object, scope = NULL, sigma.scope = NULL, nu.scope = NULL, 
              tau.scope = NULL, mu.try = TRUE, sigma.try = TRUE, 
              nu.try = TRUE, tau.try = TRUE, ...)


stepGAICAll.B(object, scope, direction = c("both", "backward", "forward"), 
               trace = T, keep = NULL, steps = 1000, scale = 0, k = 2, ...)
               
stepTGD(object, scope, newdata, direction = c("both", "backward", "forward"), 
              trace = T, keep = NULL, steps = 1000, 
              what = c("mu", "sigma", "nu", "tau"), ...)

}



\arguments{
  \item{object}{an gamlss object. This
          is used as the initial model in the stepwise search. }
  \item{scope}{defines the range of models examined in the stepwise search.
          For the function   \code{stepAIC()} this should be either a single formula, 
          or a list containing  components \code{upper} and \code{lower}, both formulae.  
          See the details for how to specify the formulae and how they are
          used.
          For the function \code{stepGAIC} the scope defines the range of models examined in the step-wise search.
          It is a list of formulas, with each formula corresponding to a term in the model. 
          A 1 in the formula allows the additional option of leaving the term out of the model entirely. +          
          }
  \item{direction}{the mode of stepwise search, can be one of \code{both},
          \code{backward}, or \code{forward}, with a default of \code{both}. If
          the \code{scope} argument is missing the default for \code{direction}
          is \code{backward}}.
  \item{trace}{if positive, information is printed during the running of
          \code{stepAIC}. Larger values may give more information on the
          fitting process.}
  \item{keep}{a filter function whose input is a fitted model object and
          the associated 'AIC' statistic, and whose output is
          arbitrary. Typically 'keep' will select a subset of the
          components of the object and return them. The default is not
          to keep anything.}
  \item{steps}{ the maximum number of steps to be considered.  The default is
          1000 (essentially as many as required).  It is typically used
          to stop the process early. }
  \item{scale}{scale is nor used in gamlss}
  \item{what}{which distribution parameter is required, default \code{what="mu"} }
  \item{k}{ the multiple of the number of degrees of freedom used for the
          penalty. Only 'k = 2' gives the genuine AIC: 'k = log(n)' is
          sometimes referred to as BIC or SBC.}
  \item{additive}{if \code{additive=TRUE} then \code{stepGAIC.CH} is used else \code{stepGAIC.CH}, default value is FALSE}
   \item{sigma.scope}{scope for \code{sigma} if different to \code{scope} in \code{stepGAICAll.A()}}
  \item{nu.scope}{scope for \code{nu} if different to \code{scope} in \code{stepGAICAll.A()}}
  \item{tau.scope}{scope for \code{tau} if different to \code{scope} in \code{stepGAICAll.A()}}
  \item{mu.try}{The default value is is TRUE, set to FALSE if no model for \code{mu} is needed}
  \item{sigma.try}{The default value is TRUE, set to FALSE if no model for \code{sigma} is needed}
  \item{nu.try}{The default value is TRUE, set to FALSE if no model for \code{nu} is needed}
  \item{tau.try}{The default value is TRUE, set to FALSE if no model for \code{tau} is needed}
   \item{newdata}{The new data set where the Test Global Deviance (TGD) will be evaluated}
  \item{\dots}{ any additional arguments to 'extractAIC'. (None are currently
          used.)  }
}
\details{
The set of models searched is determined by the \code{scope} argument.

For the function \code{stepGAIC.VR()} the right-hand-side of its \code{lower} 
component is always included in  the model, and right-hand-side of the model is included in the \code{upper} 
component.  If \code{scope} is a single formula, it specifies  the \code{upper} component,
and the \code{lower} model is empty.  If \code{scope} is missing, the initial model 
is used as the \code{upper} model.

Models specified by \code{scope} can be templates to update \code{object} as
used by \code{update.formula}.

For the function \code{stepGAIC.CH()} each of the formulas in scope specifies a 
"regimen" of candidate forms in which the particular term may enter the model. 
For example, a term formula might be 

~ x1 + log(x1) + cs(x1, df=3)

This means that x1 could either appear linearly, linearly in its logarithm, or as a smooth function estimated non-parametrically.
Every term in the model is described by such a term formula, and the final model is built up by selecting a component from each formula. 

The function \code{gamlss.scope} similar to the S \code{gam.scope()} in Chambers and Hastie (1991) can be used to create automatically
term formulae from specified data or model frames.

The supplied model object is used as the starting model, and hence there is the requirement 
that one term from each of the term formulas of the parameters be present in the formula of the distribution parameter. 
This also implies that any terms in formula of the distribution parameter not contained in any of the term formulas 
will be forced to be present in every model considered.

When the smoother used in \code{gamlss} modeling belongs to the new generation of smoothers allowing the determination of the  smoothing parameters
automatically (i.e. \code{pb()}, \code{cy()}) then the function   \code{stepGAIC.VR()} can be used for model selection (see example below).

The function \code{stepTGD} is a clone function of \code{stepGAIC.VR()} where the selection criterion is not anymore a GAIC but the Test
Global Deviance, that is, the deviance (-2log(Likelihood)) evaluated at the test data  rather than the training ones.  

}
\value{
   the stepwise-selected model is returned, with up to two additional
     components.  There is an '"anova"' component corresponding to the
     steps taken in the search, as well as a '"keep"' component if the
     'keep=' argument was supplied in the call. The '"Resid. Dev"'
     column of the analysis of deviance table refers to a constant
     minus twice the maximized log likelihood
     
   The function \code{stepGAICAll.A()} returns with a component "anovaAll" containing all the different anova tables used in the process.
}
\references{
Chambers, J. M. and Hastie, T. J. (1991). \emph{Statistical Models in S}, Chapman and Hall, London. 

Rigby, R. A. and  Stasinopoulos D. M. (2005). Generalized additive models for location, scale and shape,(with discussion), 
\emph{Appl. Statist.}, \bold{54}, part 3, pp 507-554.

Stasinopoulos D. M., Rigby R.A. and Akantziliotou C. (2006) Instructions on how to use the GAMLSS package in R.
Accompanying documentation in the current GAMLSS  help files, (see also  \url{http://www.gamlss.org/}). 

Stasinopoulos D. M. Rigby R.A. (2007) Generalized additive models for location scale and shape (GAMLSS) in R.
\emph{Journal of Statistical Software}, Vol. \bold{23}, Issue 7, Dec 2007, \url{http://www.jstatsoft.org/v23/i07}.

Venables, W. N. and Ripley, B. D. (2002) \emph{Modern Applied
     Statistics with S}. Fourth edition.  Springer.

}
\author{Mikis Stasinopoulos based on functions in MASS library and in Statistical Models in S}


\seealso{ \code{\link{gamlss.scope}} }
\examples{
\dontrun{
data(usair)
# Note default of additive=FALSE
# fitting all variables linearly 
mod1<-gamlss(y~., data=usair, family=GA)
# find the best subset for the mu
mod2<-stepGAIC(mod1)
mod2$anova
# find the best subset for sigma
mod3<-stepGAIC(mod2, what="sigma", scope=~x1+x2+x3+x4+x5+x6)
mod3$anova
# find the best model using pb() smoother 
#only three variables are used here for simplicity
mod10 <-gamlss(y~1, data=usair, family=GA)
mod20<-stepGAIC(mod10, scope=list(lower=~1, upper=~pb(x1)+pb(x2)+pb(x5)))
edf(mod20)
# x1 and x2 enter linearly
# now use the stepGAIC.CH function
# creating a scope from the usair model frame 
gs<-gamlss.scope(model.frame(y~x1+x2+x3+x4+x5+x6, data=usair))
gs 
mod4<-gamlss(y~1, data=usair, family=GA)
mod5<-stepGAIC(mod4,gs, additive=TRUE)
mod5$anova
mod6<-stepGAIC(mod5, what="sigma", scope=~x1+x2+x3+x4+x5+x6)
mod6$anova
mod6
# now stepGAICAll.A    
mod7<-stepGAICAll.A(mod4, scope=list(lower=~1,upper=~x1+x2+x3+x4+x5+x6)) 
# now  stepGAICAll.B
mod8<-stepGAICAll.B(mod4, scope=list(lower=~1,upper=~x1+x2+x3+x4+x5+x6))
# now stepTGD()
data(aep)
# sampling from the data
rand <- sample(2, dim(aep)[1], replace=TRUE, prob=c(0.6,0.4))
# the proportions in the sample
table(rand)/dim(aep)[1]
oldaep<-aep[rand==1,]
newaep<-aep[rand==2,] 
m0<-gamlss(y~ward+year+loglos, data=oldaep, family=BB)
#checking mu
m1 <- stepTGD(m0, newdata=newaep) 
#checking sigma 
m2 <- stepTGD(m0,scope=~ward+year, newdata=newaep, what="sigma") 
}
}
\keyword{regression}
