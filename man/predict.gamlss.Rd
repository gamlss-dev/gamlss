\name{predict.gamlss}
\alias{predict.gamlss}
\alias{predictAll}



\title{Extract Predictor Values and Standard Errors For New Data In a GAMLSS Model }
\description{
 \code{predict.gamlss} is the GAMLSS specific method which produce predictors for a new data set   
for a specified parameter from a GAMLSS objects.  
The \code{predict.gamlss} can be used to extract the linear predictors, fitted values and specific  terms in the model at new 
data values in the same way that the \code{predict.lm()} and \code{predict.glm()} functions can be used for 
\code{lm} or \code{glm} objects. Note that linear predictors, fitted values and specific terms in the model at the current  
data values can also be extracted using the function \code{lpred()} (which is called from predict if new data is NULL).
}
\usage{
\method{predict}{gamlss}(object, what = c("mu", "sigma", "nu", "tau"), 
               newdata = NULL, type = c("link", "response", "terms"), 
               terms = NULL, se.fit = FALSE, data = NULL, ...)
predictAll(object, newdata = NULL, type = c("response", "link", "terms"), 
               terms = NULL, se.fit = FALSE, ...)

}

\arguments{
  \item{object}{a GAMLSS fitted model}
  \item{what}{ which distribution parameter is required, default \code{what="mu"}}
  \item{newdata}{a data frame containing new values for the explanatory variables used in the model}
  \item{type}{the default, gets the linear predictor for the specified distribution parameter. 
             \code{type="response"} gets the fitted values for the parameter while \code{type="terms"} gets the fitted terms contribution}
  \item{terms}{if \code{type="terms"}, which terms to be selected (default is all terms)}
   \item{se.fit}{if TRUE the approximate standard errors of the appropriate type are extracted if exist}
  \item{data}{ the data frame used in the original fit if is not defined in the call}
  \item{\dots}{ for extra arguments }
}

\details{The predict function  assumes that the object given in \code{newdata} is a data frame containing the right x-variables
used in the model. This could possible cause problems if transformed variables are used in the fitting of the original model.
For example, let us assume that a transformation of age is needed in the model i.e. \code{nage<-age^.5}. This could be fitted as 
\code{mod<-gamlss(y~cs(age^.5),data=mydata)} or as  \code{nage<-age^.5;  mod<-gamlss(y~cs(nage), data=mydata)}.
The later could more efficient if the data are in thousands rather in hundreds. In the first case,   
the code \code{predict(mod,newdata=data.frame(age=c(34,56)))} would produce the right results.  
In the second case a new data frame has to be created containing the old data plus any new transform data. This data frame has to 
be declared in the \code{data} option. The option \code{newdata} should 
contain a data.frame with the new names and the transformed values in which prediction is required, (see the last example). }

\value{
  A vector or a matrix  depending on the options. 
}
\references{ Rigby, R. A. and  Stasinopoulos D. M. (2005). Generalized additive models for location, scale and shape,(with discussion), 
\emph{Appl. Statist.}, \bold{54}, part 3, pp 507-554.

Stasinopoulos D. M., Rigby R.A. and Akantziliotou C. (2006) Instructions on how to use the GAMLSS package in R.
Accompanying documentation in the current GAMLSS  help files, (see also  \url{http://www.gamlss.org/}). 

Stasinopoulos D. M. Rigby R.A. (2007) Generalized additive models for location scale and shape (GAMLSS) in R.
\emph{Journal of Statistical Software}, Vol. \bold{23}, Issue 7, Dec 2007, \url{http://www.jstatsoft.org/v23/i07}.
}

\author{Mikis Stasinopoulos }
\note{ This function is under development}

 
\seealso{  \code{\link{lp}},  \code{\link{lpred}}  }
\examples{
data(aids)
a<-gamlss(y~poly(x,3)+qrt, family=PO, data=aids) # 
newaids<-data.frame(x=c(45,46,47), qrt=c(2,3,4))
ap <- predict(a, newdata=newaids, type = "response")
ap
# now getting all the parameters
predictAll(a, newdata=newaids)
rm(a, ap)
data(abdom)
# transform x 
aa<-gamlss(y~cs(x^.5),data=abdom)
# predict at old values
predict(aa)[610]
# predict at new values 
predict(aa,newdata=data.frame(x=42.43))
# now transform x first 
nx<-abdom$x^.5
aaa<-gamlss(y~cs(nx),data=abdom)
# create a new data frame 
newd<-data.frame( abdom, nx=abdom$x^0.5)
# predict at old values
predict(aaa)[610]
# predict at new values 
predict(aaa,newdata=data.frame(nx=42.43^.5), data=newd)
}
\keyword{regression}% 
