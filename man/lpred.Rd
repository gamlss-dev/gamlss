\name{lpred}
\alias{lpred}
\alias{lp}

\title{Extract Linear Predictor Values and Standard Errors For A GAMLSS Model }
\description{
The function \code{lpred()} is the GAMLSS specific method which extracts the linear predictor and its (approximate) standard errors 
for a specified model parameter from a GAMLSS objects.  
The \code{lpred()} can be used to extract the predictor fitted values (and its approximate standard errors) or the contribution of specific terms in the model  
(with their approximate standard errors) in the same way that the \code{predict.lm()} and \code{predict.glm()} functions can be used for 
\code{lm} or \code{glm} objects. 
Note that \code{lpred()} extract information for the predictors of \code{mu},\code{sigma}, \code{nu} and \code{tau} at the training data values.  If predictions are  required for new data then use the 
functions \code{predict.gamlss()} or \code{predictAll()}.

The function \code{lp} extract only the linear predictor at the training data values. 
}
\usage{
lpred(obj, what = c("mu", "sigma", "nu", "tau"), parameter= NULL,
           type = c("link", "response", "terms"), 
           terms = NULL, se.fit = FALSE, ...)
lp(obj, what = c("mu", "sigma", "nu", "tau"), parameter= NULL, ... ) 
}

\arguments{
 \item{obj}{a GAMLSS fitted model}
 \item{what}{which distribution parameter is required, default \code{what="mu"} }
 \item{parameter}{equivalent to \code{what}}
 \item{type}{\code{type="link"} (the default) gets the linear predictor for the specified distribution parameter. 
             \code{type="response"} gets the fitted values for the parameter while \code{type="terms"} gets the fitted terms contribution}
 \item{terms}{if \code{type="terms"}, which terms to be selected (default is all terms) }
 \item{se.fit}{if TRUE the approximate standard errors of the appropriate type are extracted}
 \item{\dots}{for extra arguments}
}

\value{
If \code{se.fit=FALSE} a vector (or a matrix) of the appropriate \code{type} is extracted from the GAMLSS object for the given parameter in \code{what}.
If \code{se.fit=TRUE} a list containing  the appropriate \code{type}, \code{fit}, and  its (approximate) standard errors, \code{se.fit}. 
}

\references{Rigby, R. A. and  Stasinopoulos D. M. (2005). Generalized additive models for location, scale and shape,(with discussion), 
\emph{Appl. Statist.}, \bold{54}, part 3, pp 507-554.

Rigby, R. A., Stasinopoulos, D. M.,  Heller, G. Z.,  and De Bastiani, F. (2019)
	\emph{Distributions for modeling location, scale, and shape: Using GAMLSS in R}, Chapman and Hall/CRC. An older version can be found in \url{https://www.gamlss.com/}.

Stasinopoulos D. M. Rigby R.A. (2007) Generalized additive models for location scale and shape (GAMLSS) in R.
\emph{Journal of Statistical Software}, Vol. \bold{23}, Issue 7, Dec 2007, \url{https://www.jstatsoft.org/v23/i07/}.

Stasinopoulos D. M., Rigby R.A., Heller G., Voudouris V., and De Bastiani F., (2017)
\emph{Flexible Regression and Smoothing: Using GAMLSS in R},  Chapman and Hall/CRC.  

(see also \url{https://www.gamlss.com/}).
}


\author{Mikis Stasinopoulos}


\seealso{  \code{\link{predict.gamlss}} }
\examples{
data(aids)
mod<-gamlss(y~poly(x,3)+qrt, family=PO, data=aids) # 
mod.t <- lpred(mod, type = "terms", terms= "qrt")
mod.t
mod.lp <- lp(mod)
mod.lp 
rm(mod, mod.t,mod.lp)
}

\keyword{regression}% 
