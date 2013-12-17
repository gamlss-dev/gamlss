\name{VGD}
\alias{VGD}
\alias{VGD1}
\alias{VGD2}
\alias{TGD}


\title{Validation (or test) data global deviance}
\description{
These function are useful for model selection.
The   function \code{VGD} helps to validate a GAMLSS model by randomly splitting 
the data into training  and validation  sets. 
It minimies the global deviance for the training data set and 
then uses the validation set to calculate the prediction global deviance.
The function \code{VGD1} is identical to \code{VGD} but the output is a list
rather than a values as in \code{VGD}.
The function \code{VGD2} is identical to \code{VGD1} but it takes as argument
the new data, (\code{newdata}),  rather than a factor which split 
the combined data in two as  in \code{VGD} or \code{VGD1}.
The function \code{TGD} takes a fitted \code{gamlss} object and produce the test
global deviance for the new data.    
The resulting  prediction errors can be used for selecting the
distribution  of the model, the terms in the model or  degrees of freedom 
for smoothing terms.
}
\usage{
VGD(formula = NULL, sigma.formula = ~1, nu.formula = ~1, tau.formula = ~1, 
    data = NULL, family = NO, 
    control = gamlss.control(trace = FALSE), rand = NULL, ...)
VGD1(formula = NULL, sigma.formula = ~1, nu.formula = ~1, tau.formula = ~1, 
    data = NULL, family = NO, 
    control = gamlss.control(trace = FALSE), rand = NULL, ...)
VGD2(formula = NULL, sigma.formula = ~1, nu.formula = ~1, tau.formula = ~1, 
    data = NULL, family = NO, 
    control = gamlss.control(trace = FALSE), 
    newdata = NULL, ...)
TGD(object, newdata = NULL, ...)
}
\arguments{
  \item{object}{ a gamlss object}
  \item{formula}{ a gamlss formula for \code{mu} (including the response on the left)}
  \item{sigma.formula}{a formula for  \code{sigma} }
  \item{nu.formula}{ a formula for  \code{nu} }
  \item{tau.formula}{ a formula for  \code{tau} }
  \item{data}{the data set used for the fitting}
  \item{newdata}{the data set to be used for validation ot test}
  \item{family}{a \code{\link[gamlss.dist]{gamlss.family}} object }
  \item{control}{\code{gamlss.control} to be passed to \code{gamlss}}
  \item{rand}{a random vector of one and two indicating whether is the training set (1) or the 
    validation set (2) i.e. created in advance using something like 
      \code{rand <- sample(2, N, replace=T, prob=c(0.6,0.4))} where N is the length of the data}
  \item{\dots}{ for extra arguments to be passed in the \code{gamlss} fit }
}

\value{
 The function \code{VGD} returns a validated global deviance, 
 that is,  the global deviance evaluated at the validation data. 
 The functions \code{VGD1}, \code{VGD2} and \code{TGD} are returning a list
 with component \code{oldGD}, \code{newGD}, \code{oldPE}, \code{newPE}. 
 \code{oldGD} is the deviance of the training data, \code{newGD} is the deviance for 
 the validation or (test) data. \code{oldPE} is  
 the mean prediction error for the old data 
 (\code{oldGD} divided by the number of observations in the training data).
 \code{newPE}is the mean prediction error for the new data, 
 (\code{newGD} divided by the number of observations in the validation (or test) data). 
}
\references{Rigby, R. A. and  Stasinopoulos D. M. (2005). Generalized additive models for location, scale and shape,(with discussion), 
\emph{Appl. Statist.}, \bold{54}, part 3, pp 507-554.

Stasinopoulos D. M., Rigby R.A. and Akantziliotou C. (2006) Instructions on how to use the GAMLSS package in R.
Accompanying documentation in the current GAMLSS  help files, (see also  \url{http://www.gamlss.org/}).

Stasinopoulos D. M. Rigby R.A. (2007) Generalized additive models for location scale and shape (GAMLSS) in R.
\emph{Journal of Statistical Software}, Vol. \bold{23}, Issue 7, Dec 2007, \url{http://www.jstatsoft.org/v23/i07}.
}
\author{ Mikis Stasinopoulos \email{d.stasinopoulos@londonmet.ac.uk} and Bob Rigby \email{r.rigby@londonmet.ac.uk}  }

\seealso{\code{\link[gamlss.dist]{gamlss.family}}, \code{\link{gamlss}},  \code{\link{deviance.gamlss}}}
\examples{
data(abdom)
# generate the random split of the data
rand <- sample(2, 610, replace=TRUE, prob=c(0.6,0.4))
# the proportions in the sample
table(rand)/610
# VGD needs a factor
VGD(y~pb(x,df=2),sigma.formula=~pb(x,df=1), data=abdom, family=LO, rand=rand)
# VGD1 needs a factor and has different output
VGD1(y~pb(x,df=2),sigma.formula=~pb(x,df=1), data=abdom, family=LO, rand=rand)
# using VGD2 needs two different sets of data 
VGD2(y~pb(x,df=2),sigma.formula=~pb(x,df=1), data=abdom[rand==1,], family=LO,
      newdata=abdom[rand==2,] )
#fit in the old data
olddata<-abdom[rand==1,]
m1<-gamlss(y~pb(x,df=2),sigma.formula=~pb(x,df=1), data=olddata, family=LO)
# get the global deviance in the new data
newdata<-abdom[rand==2,] 
TGD(m1, newdata=newdata)
}
\keyword{regression}% 
