\name{lpred}
\alias{lpred}
\alias{lp}

\title{Extract Linear Predictor Values and Standard Errors For A GAMLSS Model }
\description{
\code{lpred} is the GAMLSS specific method which extracts the linear predictor and its (approximate) standard errors 
for a specified parameter from a GAMLSS objects.  
The \code{lpred} can be also used to extract the fitted values (with its approximate standard errors) or specific  terms in the model  
(with its approximate standard errors) in the same way that the \code{predict.lm()} and \code{predict.glm()} functions can be used for 
\code{lm} or \code{glm} objects. 
The function \code{lp} extract only the linear predictor. If prediction is required for new data values then use the 
function \code{predict.gamlss()}.
}
\usage{
lpred(obj, what = c("mu", "sigma", "nu", "tau"), 
           type = c("link", "response", "terms"), 
           terms = NULL, se.fit = FALSE, ...)
lp(obj, what = "mu", ...)
}

\arguments{
 \item{obj}{a GAMLSS fitted model}
 \item{what}{which distribution parameter is required, default \code{what="mu"} }
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

Stasinopoulos D. M., Rigby R.A. and Akantziliotou C. (2006) Instructions on how to use the GAMLSS package in R.
Accompanying documentation in the current GAMLSS  help files, (see also  \url{http://www.gamlss.org/}). 

Stasinopoulos D. M. Rigby R.A. (2007) Generalized additive models for location scale and shape (GAMLSS) in R.
\emph{Journal of Statistical Software}, Vol. \bold{23}, Issue 7, Dec 2007, \url{http://www.jstatsoft.org/v23/i07}.
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
