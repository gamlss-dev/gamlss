\name{stepGAIC}

\alias{stepGAIC
}
\alias{stepGAIC.CH}
\alias{stepGAIC.VR}
\alias{stepGAICAll.A}
\alias{stepGAICAll.B}
\alias{drop1All}
\alias{add1All}

\title{Choose a model by GAIC in a Stepwise Algorithm}

\description{

The function \code{stepGAIC()} performs stepwise model selection using a Generalized Akaike Information Criterion (GAIC). It is based on the function \code{stepAIC()} given in the library MASS of Venables and Ripley (2002). The function has been changed recently to allow parallel computation.  The parallel computations are similar to the ones performed in the function \code{boot()} of the \pkg{boot} package. Note that since version 4.3-5 of \pkg{gamlss} the  \code{stepGAIC()} do not have the option of using the function \code{stepGAIC.CH} through the argument \code{additive}. 

Note that \code{stepGAIC()}  is relying to the \code{dropterm()} and \code{addterm()} methods applied to gamlss objects.    \code{drop1()} and \code{add1()} are equivalent methods  to the  \code{dropterm()} and \code{addterm()} respectively but with different default arguments (see the examples). 

The function \code{stepGAIC.VR()} is the old version  of \code{stepGAIC()} with no parallel computations.    

The function  \code{stepGAIC.CH} is based on the S function \code{step.gam()} (see Chambers and Hastie (1991)) and it is more suited for model with smoothing additive terms when the degrees of freedom for smoothing are fixed in advance. This is something which rarely used these days, as most of the smoothing functions allow the calculations of the smoothing parameter,  see for example the additive function \code{pb()}).  

The functions \code{stepGAIC.VR()} and \code{stepGAIC.CH()} have been adapted to work with  gamlss objects and  the main difference is the \code{scope} argument, see below. 


While the functions \code{stepGAIC()} is used to build models for individual parameters of the distribution  of the response variable, the functions \code{stepGAICAll.A()} and \code{stepGAICAll.A()}  are building  models for all  the parameters.

The functions  \code{stepGAICAll.A()} and \code{stepGAICAll.B()} are based on the \code{stepGAIC()} function but use  different  strategies for selecting a appropriate final model.    

\code{stepGAICAll.A()} has the following strategy: 
 
Strategy A:
 
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


\code{stepGAICAll.B()} uses the same procedure as the function  \code{stepGAIC()} but each term in the scope is fitted  to all the parameters of the distribution, rather than the one specified  by the argument \code{what} of \code{stepGAIC()}. 
The \code{stepGAICAll.B()} relies on the \code{add1All()} and \code{drop1All()} functions for the selection of variables.
 

}
\usage{

stepGAIC(object, scope, direction = c("both", "backward", "forward"), 
          trace = TRUE, keep = NULL, steps = 1000, scale = 0, 
          what = c("mu", "sigma", "nu", "tau"), parameter= NULL, k = 2, 
          parallel = c("no", "multicore", "snow"), ncpus = 1L, cl = NULL, 
          ...)

stepGAIC.VR(object, scope, direction = c("both", "backward", "forward"), 
         trace = TRUE, keep = NULL, steps = 1000, scale = 0, 
         what = c("mu", "sigma", "nu", "tau"), parameter= NULL, k = 2, 
         ...)

stepGAIC.CH(object, scope = gamlss.scope(model.frame(object)), 
            direction = c("both", "backward", "forward"), trace = TRUE, 
            keep = NULL, steps = 1000, what = c("mu", "sigma", "nu", "tau"),
            parameter= NULL, k = 2, ...)

stepGAICAll.A(object, scope = NULL, sigma.scope = NULL, nu.scope = NULL, 
              tau.scope = NULL, mu.try = TRUE, sigma.try = TRUE, 
              nu.try = TRUE, tau.try = TRUE,  direction = NULL,  
              parallel = c("no", "multicore", "snow"), ncpus = 1L, 
              cl = NULL,  ...)

stepGAICAll.B(object, scope, direction = c("both", "backward", "forward"), 
              trace = T, keep = NULL, steps = 1000, scale = 0, k = 2, 
              parallel = c("no", "multicore", "snow"), ncpus = 1L, 
              cl = NULL, ...) 
               
drop1All(object, scope, test = c("Chisq", "none"), k = 2, sorted = FALSE, 
              trace = FALSE, parallel = c("no", "multicore", "snow"), 
              ncpus = 1L, cl = NULL, ...)
              
add1All(object, scope, test = c("Chisq", "none"), k = 2, sorted = FALSE, 
              trace = FALSE, parallel = c("no", "multicore", "snow"), 
              ncpus = 1L, cl = NULL, ...)              
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
  \item{parameter}{equivalent to \code{what}}
  \item{k}{ the multiple of the number of degrees of freedom used for the
          penalty. Only 'k = 2' gives the genuine AIC: 'k = log(n)' is
          sometimes referred to as BIC or SBC.}
   \item{parallel}{The type of parallel operation to be used (if any). If missing, the default is "no".}
   \item{ncpus}{integer: number of processes to be used in parallel operation: typically     one would chose this to the number of available CPUs.}
  \item{cl}{An optional parallel or snow cluster for use if \code{parallel = "snow"}. If not supplied, a cluster on the local machine is created for the duration of the call.
}
   \item{sigma.scope}{scope for \code{sigma} if different to \code{scope} in \code{stepGAICAll.A()}}
  \item{nu.scope}{scope for \code{nu} if different to \code{scope} in \code{stepGAICAll.A()}}
  \item{tau.scope}{scope for \code{tau} if different to \code{scope} in \code{stepGAICAll.A()}}
  \item{mu.try}{The default value is is TRUE, set to FALSE if no model for \code{mu} is needed}
  \item{sigma.try}{The default value is TRUE, set to FALSE if no model for \code{sigma} is needed}
  \item{nu.try}{The default value is TRUE, set to FALSE if no model for \code{nu} is needed}
  \item{tau.try}{The default value is TRUE, set to FALSE if no model for \code{tau} is needed}
  \item{test}{whether to print the chi-square test or not}
  \item{sorted}{whether to sort the results}        
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

When the smoother used in \code{gamlss} modelling belongs to the new generation of smoothers allowing the determination of the  smoothing parameters
automatically (i.e. \code{pb()}, \code{cy()}) then the function   \code{stepGAIC.VR()} can be used for model selection (see example below).

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


Rigby, R. A., Stasinopoulos, D. M.,  Heller, G. Z.,  and De Bastiani, F. (2019)
	\emph{Distributions for modeling location, scale, and shape: Using GAMLSS in R}, Chapman and Hall/CRC. An older version can be found in \url{https://www.gamlss.com/}.

Stasinopoulos D. M. Rigby R.A. (2007) Generalized additive models for location scale and shape (GAMLSS) in R.
\emph{Journal of Statistical Software}, Vol. \bold{23}, Issue 7, Dec 2007, \url{https://www.jstatsoft.org/v23/i07/}.

Stasinopoulos D. M., Rigby R.A., Heller G., Voudouris V., and De Bastiani F., (2017)
\emph{Flexible Regression and Smoothing: Using GAMLSS in R},  Chapman and Hall/CRC.  

(see also \url{https://www.gamlss.com/}).

Venables, W. N. and Ripley, B. D. (2002) \emph{Modern Applied
     Statistics with S}. Fourth edition.  Springer.

}
\author{Mikis Stasinopoulos based on functions in MASS library and in Statistical Models in S}


\seealso{ \code{\link{gamlss.scope}} }
\examples{
\dontrun{
data(usair)
# -----------------------------------------------------------------------------
#  null model 
mod0<-gamlss(y~1, data=usair, family=GA)
#  all the explanatotory variables x1:x6 fitted linearly 
mod1<-gamlss(y~., data=usair, family=GA)
#-------------------------------------------------------------------------------
# droping terms 
dropterm(mod1)
# with chi-square information
drop1(mod1)
# for parallel computations use something like 
nC <- detectCores()
drop1(mod1, parallel="snow",  ncpus=nC)
drop1(mod1, parallel="multicore",  ncpus=nC)
#------------------------------------------------------------------------------
# adding terms
addterm(mod0, scope=as.formula(paste("~", paste(names(usair[-1]),
                  collapse="+"),sep="")))
# with chi-square information
add1(mod0, scope=as.formula(paste("~", paste(names(usair[-1]),
                  collapse="+"),sep="")))
# for parallel computations 
nC <- detectCores()
add1(mod0, scope=as.formula(paste("~", paste(names(usair[-1]),
                  collapse="+"),sep="")), parallel="snow",  ncpus=nC)
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
# stepGAIC 
# find the best subset for the mu
mod2 <- stepGAIC(mod1)
mod2$anova
#--------------------------------------------------------------
# for parallel computations 
mod21 <- stepGAIC(mod1, , parallel="snow",  ncpus=nC)
#--------------------------------------------------------------
# find the best subset for sigma
mod3<-stepGAIC(mod2, what="sigma", scope=~x1+x2+x3+x4+x5+x6)
mod3$anova
#--------------------------------------------------------------
# find the best model using pb() smoother 
#only three variables are used here for simplicity
mod20<-stepGAIC(mod0, scope=list(lower=~1, upper=~pb(x1)+pb(x2)+pb(x5)))
edf(mod20)
# note that x1 and x2 enter linearly
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
# the stepGAIC.CH function (no parallel here)
# creating a scope from the usair model frame 
gs<-gamlss.scope(model.frame(y~x1+x2+x3+x4+x5+x6, data=usair))
gs 
mod5<-stepGAIC.CH(mod0,gs)
mod5$anova
#-----------------------------------------------------------------------------=-
#------------------------------------------------------------------------------
# now stepGAICAll.A    
mod7<-stepGAICAll.A(mod0, scope=list(lower=~1,upper=~x1+x2+x3+x4+x5+x6)) 
#-----------------------------------------------------------------------------=-
#------------------------------------------------------------------------------
# now  stepGAICAll.B
drop1All(mod1, parallel="snow",  ncpus=nC)
add1All(mod0, scope=as.formula(paste("~", paste(names(usair[-1]),
                  collapse="+"))), parallel="snow",  ncpus=nC)
mod8<-stepGAICAll.B(mod0, scope=list(lower=~1,upper=~x1+x2+x3+x4+x5+x6))
#-----------------------------------------------------------------------------=-
#------------------------------------------------------------------------------
}
}
\keyword{regression}
