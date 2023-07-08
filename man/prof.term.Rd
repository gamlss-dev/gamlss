\name{prof.term}
\alias{prof.term}

\title{Plotting the Profile: deviance or information criterion for one of the terms (or hyper-parameters) in a GAMLSS model }
\description{
This functions plots the profile deviance for a chosen parameter included in the linear predictor of any of the \code{mu},\code{sigma}, \code{nu} or \code{tau} models so profile confidence intervals can be obtained. In can also be used to plot the profile of a specified information criterion for any hyper-parameter when smooth additive terms are used.  
}
\usage{
prof.term(model = NULL, criterion = c("GD", "GAIC"), penalty = 2.5, 
          other = NULL, min = NULL, max = NULL, step = NULL, 
          length = 7, xlabel = NULL, plot = TRUE, perc = 95, 
          start.prev = TRUE, col="darkgreen")
}

\arguments{
  \item{model}{this is a GAMLSS model, e.g. \cr
                     \code{model=gamlss(y~cs(x,df=this),sigma.fo=~cs(x,df=3),data=abdom)}, where \code{this} indicates the (hyper)parameter to be profiled}
  \item{criterion}{whether global deviance ("GD") or information criterion ("GAIC") is profiled. The default is global deviance \code{criterion="GD"}}
  \item{penalty}{ The penalty value if information criterion is used in \code{criterion}, default \code{penalty=2.5}  }
  \item{other}{this can be used to evaluate an expression before the actual fitting of the model (Make sure that those expressions are well define in the global environment)}
  \item{min}{the minimum value for the parameter e.g. \code{min=1} }
  \item{max}{the maximum value for the parameter e.g.  \code{max=20}}
  \item{step}{how often to evaluate the global deviance (defines the step length of the grid for the parameter) e.g. \code{step=1}}
  \item{length}{if the step is left NULL then \code{length} is considered for evaluating the grid for the parameter. It has a default value of 11}
  \item{xlabel}{if a label for the axis is required }
  \item{plot}{whether to plot, \code{plot=TRUE} the resulting profile deviance (or GAIC)}
  \item{perc}{what \% confidence interval is required}
  \item{start.prev}{whether to start from the previous fitted model parameters values or not (default is TRUE)}
  \item{col}{the color of the profile line}
}
\details{
 This function can be use to provide likelihood based confidence intervals for a parameter involved in terms in the linear predictor(s). 
 These confidence intervals are more accurate than the ones obtained from the parameters' standard errors. 
 The function can also be used to plot a profile information criterion (with a given penalty) against a hyper-parameter. This can be used to check the uniqueness in hyper-parameter determination using for example \code{find.df}. 
}
\value{
  Return a profile plot (if the argument \code{plot=TRUE}) and an \code{ProfLikelihood.gamlss} object if saved.  The object contains:
  
  \item{values}{the values at the grid where the parameter was evaluated}
  \item{fun}{the function which approximates the points using splines}
   \item{min}{the minimum values in the grid}
    \item{max}{the maximum values in the grid}
     \item{max.value}{the value of the parameter maximising the Profile deviance (or GAIC) }
      \item{CI}{the profile confidence interval (if global deviance is used)}
       \item{criterion}{which criterion was used}
}
\references{
Rigby, R. A. and  Stasinopoulos D. M. (2005). Generalized additive models for location, scale and shape,(with discussion), 
\emph{Appl. Statist.}, \bold{54}, part 3, pp 507-554.

Rigby, R. A., Stasinopoulos, D. M.,  Heller, G. Z.,  and De Bastiani, F. (2019)
	\emph{Distributions for modeling location, scale, and shape: Using GAMLSS in R}, Chapman and Hall/CRC. An older version can be found in \url{https://www.gamlss.com/}.

Stasinopoulos D. M. Rigby R.A. (2007) Generalized additive models for location scale and shape (GAMLSS) in R.
\emph{Journal of Statistical Software}, Vol. \bold{23}, Issue 7, Dec 2007, \url{https://www.jstatsoft.org/v23/i07/}.

Stasinopoulos D. M., Rigby R.A., Heller G., Voudouris V., and De Bastiani F., (2017)
\emph{Flexible Regression and Smoothing: Using GAMLSS in R},  Chapman and Hall/CRC.  

(see also \url{https://www.gamlss.com/}).
}
\author{ Mikis Stasinopoulos \email{d.stasinopoulos@londonmet.ac.uk} and Bob Rigby}
\section{Warning }{A dense grid (i.e. small step) evaluation of the global deviance can take a long time, so start with a sparse grid (i.e. large step) 
and decrease  gradually the step length for more accuracy. }  
\seealso{\code{\link{gamlss}},  \code{\link{prof.dev}} }
\examples{
data(aids)
# fitting a linear model
gamlss(y~x+qrt,family=NBI,data=aids)
# testing the linear beta parameter
mod<-quote(gamlss(y ~ offset(this * x) + qrt, data = aids, family = NBI))
prof.term(mod, min=0.06, max=0.11)
# find the hyper parameter using cubic splines smoothing
mod1<-quote(gamlss(y ~ cs(x,df=this) + qrt, data = aids, family = NBI))
prof.term(mod1, min=1, max=15, step=1, criterion="GAIC", penalty=log(45))
# find a break point in x
mod2 <- quote(gamlss(y ~ x+I((x>this)*(x-this))+qrt,family=NBI,data=aids))
prof.term(mod2, min=1, max=45, step=1, criterion="GD")
rm(mod,mod1,mod2)
}
\keyword{regression}% 
