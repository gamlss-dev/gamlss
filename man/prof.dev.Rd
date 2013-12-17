\name{prof.dev}
\alias{prof.dev}

\title{Plotting the Profile Deviance for one of the Parameters in a GAMLSS model}
\description{
This functions plots the profile deviance of one of the (four) parameters in a GAMLSS model. It can be used if one 
of the parameters \code{mu}, \code{sigma}, \code{nu} or \code{tau} is a constant (not a function of explanatory variables) to obtain 
a profile confidence intervals. 
}
\usage{
prof.dev(object, which = NULL, min = NULL, max = NULL, 
        step = NULL, length = 7, startlastfit = TRUE, 
        plot = TRUE, perc = 95, ...)
}

\arguments{
  \item{object}{A fitted GAMLSS model}
  \item{which}{which parameter to  get the profile deviance e.g. \code{which="tau"}}
  \item{min}{the minimum value for the parameter e.g. \code{min=1}}
  \item{max}{the maximum value for the parameter e.g.  \code{max=20}}
  \item{step}{how often to evaluate the global deviance (defines the step length of the grid for the parameter) e.g. \code{step=1} }
  \item{length}{the length if step is not set, default equal 7}
  \item{startlastfit}{whether to start fitting from the last fit or not, default value is \code{startlastfit=TRUE}  }
   \item{plot}{whether to plot, \code{plot=TRUE} or save the results, \code{plot=FALSE}  }
  \item{perc}{what \% confidence interval is required }
  \item{\dots}{for extra arguments}
}
\details{
This function can be use to provide likelihood based confidence intervals for a parameter for which a constant model (i.e. no explanatory model) is fitted and
consequently for checking the adequacy of a particular values of the parameter. This can be used to check the adequacy of one distribution (e.g. Box-Cox Cole and Green) 
nested within another (e.g. Box-Cox power exponential). For example one can test whether a Box-Cox Cole and Green (Box-Cox-normal) distribution 
or a Box-Cox power exponential is appropriate by plotting the profile of the parameter \code{tau}. 
A profile deviance showing support for  \code{tau=2} indicates adequacy of the Box-Cox Cole and Green (i.e. Box-Cox normal) distribution.      
}
\value{
 Return a profile plot (if the argument \code{plot=TRUE}) and an \code{ProfLikelihood.gamlss} object if saved.  The object contains:
  
  \item{values}{the values at the grid where the parameter was evaluated}
  \item{fun}{the function which approximates the points using splines}
   \item{min}{the minimum values in the grid}
    \item{max}{te maximum values in the grid}
     \item{max.value}{the value of the parameter maximising the Profile deviance (or GAIC) }
      \item{CI}{the profile confidence interval (if global deviance is used)}
       \item{criterion}{which criterion was used}
}
\references{Rigby, R. A. and  Stasinopoulos D. M. (2005). Generalized additive models for location, scale and shape,(with discussion), 
\emph{Appl. Statist.}, \bold{54}, part 3, pp 507-554.

Stasinopoulos D. M., Rigby R.A. and Akantziliotou C. (2006) Instructions on how to use the GAMLSS package in R.
Accompanying documentation in the current GAMLSS  help files, (see also  \url{http://www.gamlss.org/}).  

Stasinopoulos D. M. Rigby R.A. (2007) Generalized additive models for location scale and shape (GAMLSS) in R.
\emph{Journal of Statistical Software}, Vol. \bold{23}, Issue 7, Dec 2007, \url{http://www.jstatsoft.org/v23/i07}.
}


\author{Calliope Akantziliotou, Mikis Stasinopoulos \email{d.stasinopoulos@londonmet.ac.uk} and Bob Rigby \email{r.rigby@londonmet.ac.uk} }
\section{Warning }{A dense grid (i.e. small step) evaluation of the global deviance can take a long time, so start with a sparse grid (i.e. large step) 
and decrease  gradually the step length for more accuracy. }  
\seealso{ \code{\link{gamlss}},  \code{\link{prof.term}} }
\examples{
\dontrun{
data(abdom)
h<-gamlss(y~pb(x), sigma.formula=~pb(x), family=BCT, data=abdom) 
prof.dev(h,"nu",min=-2.000,max=2)
rm(h)}
}
\keyword{regression}% 
