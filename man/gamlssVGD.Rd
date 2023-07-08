\name{gamlssVGD}
\alias{gamlssVGD}
\alias{VGD}
\alias{getTGD}
\alias{TGD}
\alias{gamlssCV}
\alias{CV}
\alias{drop1TGD}
\alias{add1TGD}
\alias{stepTGD}
\alias{stepTGDAll.A}



%- Also NEED an '\alias' for EACH other topic documented here.
\title{A  Set of Functions for selecting Models using Validation or Test Data Sets and Cross Validation 
}
\description{
This is a set of function useful for selecting appropriate models. 

The functions \code{gamlssVGD}, \code{VGD}, \code{getTGD}, \code{TGD} can be used when a subset of the data is used for validation or testing. 

The function \code{stepVGD()} is a stepwise procedure for selecting an appropriate model for any of the parameters of the model minimising the test global deviance. The function   \code{stepVGDAll.A()} can select a model using strategy A for all the parameters. 

The functions \code{gamlssCV}, \code{CV} can be used for a k-fold cross validation.  


}
\usage{
gamlssVGD(formula = NULL, sigma.formula = ~1, nu.formula = ~1, 
          tau.formula = ~1, data = NULL, family = NO, 
          control = gamlss.control(trace = FALSE), 
          rand = NULL, newdata = NULL, ...)
          
VGD(object, ...)          

getTGD(object, newdata = NULL, ...)

TGD(object, ...)  

gamlssCV(formula = NULL, sigma.formula = ~1, nu.formula = ~1, 
         tau.formula = ~1, data = NULL, family = NO, 
         control = gamlss.control(trace = FALSE), 
         K.fold = 10, set.seed = 123, rand = NULL, 
         parallel = c("no", "multicore", "snow"), 
         ncpus = 1L, cl = NULL, ...)

CV(object, ...)

drop1TGD(object, scope, newdata, parameter = c("mu", "sigma", "nu", "tau"), 
         sorted = FALSE, trace = FALSE, 
         parallel = c("no", "multicore", "snow"), 
         ncpus = 1L, cl = NULL, ...)
         
add1TGD(object, scope, newdata, parameter = c("mu", "sigma", "nu", "tau"), 
        sorted = FALSE, trace = FALSE, 
        parallel = c("no", "multicore", "snow"), 
        ncpus = 1L, cl = NULL, ...)

stepTGD(object, scope, newdata, 
        direction = c("both", "backward", "forward"),
        trace = TRUE, keep = NULL, steps = 1000, 
        parameter = c("mu", "sigma", "nu", "tau"), 
        parallel = c("no", "multicore", "snow"), 
        ncpus = 1L, cl = NULL, ...)
        
stepTGDAll.A(object, scope = NULL, newdata = NULL, 
        steps = 1000, sigma.scope = NULL, nu.scope = NULL, 
        tau.scope = NULL, mu.try = TRUE, sigma.try = TRUE, 
        nu.try = TRUE, tau.try = TRUE,
        parallel = c("no", "multicore", "snow"), 
        ncpus = 1L, cl = NULL, ...)        

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{A \code{gamlss} \code{mu} formula.}
  \item{sigma.formula}{Formula for \code{sigma}.}
  \item{nu.formula}{Formula for \code{nu}.}
  \item{tau.formula}{Formula for \code{tau}.}
  \item{data}{The data frame required for the fit.}
  \item{family}{The \code{gamlss.family} distribution.}
  \item{control}{The control for fitting the gamlss model.}
  \item{rand}{For \code{gamlssVGD} a variable with values 1 (for fitting) and 2 (for predicting). For \code{gamlssCV} a variable with k values indicating the cross validation sets.}
  \item{newdata}{The new data set (validation or test) for prediction.}
  \item{object}{A relevant R object.}
  \item{scope}{defines the range of models examined in the stepwise selection similar to \code{stepGAIC()} where you can see examples}
  \item{sigma.scope}{defines the range of models examined in the stepwise selection for \code{sigma}} 
  \item{nu.scope}{defines the range of models examined in the stepwise selection for \code{nu}} 
  \item{tau.scope}{defines the range of models examined in the stepwise selection for \code{tau}} 
  \item{mu.try}{whether should try fitting models for \code{mu}} 
  \item{sigma.try}{whether should try fitting models for \code{sigma}} 
  \item{nu.try}{whether should try fitting models for \code{nu}} 
  \item{tau.try}{whether should try fitting models for \code{tau}} 
  \item{parameter}{which distribution parameter is required, default \code{what="mu"}}
  \item{sorted}{should the results be sorted on the value of TGD}
  \item{trace}{f \code{TRUE} additional information may be given on the fits as they are tried.}
  \item{direction}{The mode of stepwise search, can be one of \code{both}, \code{backward}, or \code{forward}, with a default of \code{both}. If the scope argument is missing the default for direction is backward}
  \item{keep}{see \code{stepGAIC()} for explanation}
  \item{steps}{the maximum number of steps to be considered. The default is 1000.}
  
  \item{K.fold}{the number of subsets of the data used} 
   \item{set.seed}{the seed to be used in creating \code{rand}}
   \item{parallel}{The type of parallel operation to be used (if any). If missing, the default is "no".}
  \item{ncpus}{integer: number of processes to be used in parallel operation: typically one would chose this to the number of available CPUs.}
  \item{cl}{An optional parallel or snow cluster for use if \code{parallel = "snow"}. If not supplied, a cluster on the local machine is created for the duration of the call.
}
  \item{\dots}{further arguments to be pass in the gamlss fit}
}
\details{
The function \code{gamlssVGD()} fits a gamlss model to the training data set determined by the arguments \code{rand} or \code{newdata}. The results is a \code{gamlssVGD} objects which contains the gamlss fit to the training data plus three extra components: i) \code{VGD} the global deviance applied to the validation data sets. ii) \code{predictError} which is \code{VGD}
divided with the number of observations in the validation data set and iii) \code{residVal} the residuals for the  validation data set. 

The function \code{VGD()} extract the validated global deviance from one or more fitted  \code{gamlssVGD} objects and can be used foe model comparison. 

The function \code{getTGD()} operates different from the function \code{gamlssVGD()}. It assumes that the users already have fitted models using \code{gamlss()} and now he/she wants to evaluate the global deviance at a new (validation or test) data set. 

The function \code{TGD()}  extract the validated/test global deviance  from one or more fitted  \code{gamlssTGD} objects and can be use to compare models.

The \code{gamlssCV()} performs a k-fold cross validation on a gamlss models. 

The function \code{CV()} extract the cross validated  global deviance from one or more fitted  \code{gamlssCV} objects and can be use to compare models.

The functions \code{add1TGD()}, \code{drop1TGD()} and \code{stepTGD} behave similar to  \code{add1()}, \code{drop1()} and \code{stepGAIC()} functions respectively  but they used  validation or test deviance as the selection criterion rather than the GAIC.



}
\value{
A fitted models of a set of global deviances.
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
\author{Mikis Stasinopoulos}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{\link{stepGAIC}}}
\examples{
data(abdom)
# generate the random split of the data
rand <- sample(2, 610, replace=TRUE, prob=c(0.6,0.4))
# the proportions in the sample
table(rand)/610
olddata<-abdom[rand==1,] # training data
newdata<-abdom[rand==2,] # validation data
#------------------------------------------------------------------------------
#  gamlssVGD
#-------------------------------------------------------------------------------
# Using rand
v1 <- gamlssVGD(y~pb(x,df=2),sigma.formula=~pb(x,df=1), data=abdom, family=NO, 
              rand=rand)
v2 <- gamlssVGD(y~pb(x,df=2),sigma.formula=~pb(x,df=1), data=abdom, family=LO, 
              rand=rand)
v3 <- gamlssVGD(y~pb(x,df=2),sigma.formula=~pb(x,df=1), data=abdom, family=TF, 
              rand=rand)
VGD(v1,v2,v3)
#-------------------------------------------------------------------------------
\dontrun{
#-------------------------------------------------------------------------------
# using two data set
v11 <- gamlssVGD(y~pb(x,df=2),sigma.formula=~pb(x,df=1), data=olddata, 
                  family=NO, newdata=newdata)
v12 <- gamlssVGD(y~pb(x,df=2),sigma.formula=~pb(x,df=1), data=olddata, 
                 family=LO, newdata=newdata)
v13 <- gamlssVGD(y~pb(x,df=2),sigma.formula=~pb(x,df=1), data=olddata, 
                 family=TF, newdata=newdata)
VGD(v11,v12,v13)
#-------------------------------------------------------------------------------
# function getTGD
#-------------------------------------------------------------------------------
# fit gamlss models first
g1 <- gamlss(y~pb(x,df=2),sigma.formula=~pb(x,df=1), data=olddata, family=NO)
g2 <- gamlss(y~pb(x,df=2),sigma.formula=~pb(x,df=1), data=olddata, family=LO)
g3 <- gamlss(y~pb(x,df=2),sigma.formula=~pb(x,df=1), data=olddata, family=TF)
# and then use 
gg1 <-getTGD(g1, newdata=newdata)
gg2 <-getTGD(g2, newdata=newdata)
gg3 <-getTGD(g3, newdata=newdata)

TGD(gg1,gg2,gg3)
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# function gamlssCV
#-------------------------------------------------------------------------------
set.seed(123)
rand1 <- sample (10 , 610, replace=TRUE)
g1 <- gamlssCV(y~pb(x,df=2),sigma.formula=~pb(x,df=1), data=abdom, family=NO, 
               rand=rand1)
g2 <- gamlssCV(y~pb(x,df=2),sigma.formula=~pb(x,df=1), data=abdom, family=LO, 
               rand=rand1)
g3 <- gamlssCV(y~pb(x,df=2),sigma.formula=~pb(x,df=1), data=abdom, family=TF, 
               rand=rand1)
CV(g1,g2,g3)
CV(g1)
# using parallel 
set.seed(123)
rand1 <- sample (10 , 610, replace=TRUE)
nC <- detectCores()

system.time(g21 <- gamlssCV(y~pb(x,df=2), sigma.formula=~pb(x,df=1), data=abdom,
             family=NO, rand=rand1,parallel = "no", ncpus = nC ))

system.time(g22 <- gamlssCV(y~pb(x,df=2), sigma.formula=~pb(x,df=1), data=abdom,
             family=LO, rand=rand1,parallel = "multicore", ncpus = nC ))

system.time(g23 <- gamlssCV(y~pb(x,df=2), sigma.formula=~pb(x,df=1), data=abdom,
             family=TF, rand=rand1,parallel = "snow", ncpus = nC ))


CV(g21,g22,g23)
#-------------------------------------------------------------------------------
# functions add1TGD() drop1TGD() and stepTGD()
#-------------------------------------------------------------------------------
# the data
data(rent)
rand <- sample(2, dim(rent)[1], replace=TRUE, prob=c(0.6,0.4))
# the proportions in the sample
table(rand)/dim(rent)[1]
oldrent<-rent[rand==1,] # training set
newrent<-rent[rand==2,] # validation set

# null model
v0 <- gamlss(R~1, data=oldrent, family=GA)
# complete model
v1 <- gamlss(R~pb(Fl)+pb(A)+H+loc, sigma.fo=~pb(Fl)+pb(A)+H+loc, 
             data=oldrent, family=GA)

# drop1TGDP
system.time(v3<- drop1TGD(v1, newdata=newrent,  parallel="no"))
system.time(v4<- drop1TGD(v1, newdata=newrent,  parallel="multicore", 
                          ncpus=nC) )
system.time(v5<- drop1TGD(v1, newdata=newrent,  parallel="snow", ncpus=nC))
cbind(v3,v4,v5)

# add1TGDP
system.time(d3<- add1TGD(v0,scope=~pb(Fl)+pb(A)+H+loc, newdata=newrent,  
                       parallel="no"))
system.time(d4<- add1TGD(v0,scope=~pb(Fl)+pb(A)+H+loc, newdata=newrent,  
                        parallel="multicore", ncpus=nC) )
system.time(d5<- add1TGD(v0, scope=~pb(Fl)+pb(A)+H+loc,newdata=newrent,  
                        parallel="snow", ncpus=nC))

# stepTGD
system.time(d6<- stepTGD(v0, scope=~pb(Fl)+pb(A)+H+loc,newdata=newrent))
system.time(d7<- stepTGD(v0, scope=~pb(Fl)+pb(A)+H+loc,newdata=newrent,
                         parallel="multicore", ncpus=nC))
system.time(d8<- stepTGD(v0, scope=~pb(Fl)+pb(A)+H+loc,newdata=newrent,
                         parallel="snow", ncpus=nC))
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{regression}
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
