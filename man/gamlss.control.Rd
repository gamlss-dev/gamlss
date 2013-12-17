\name{gamlss.control}
\alias{gamlss.control}

\title{Auxiliary for Controlling GAMLSS Fitting
}
\description{
 Auxiliary function as user interface for \code{gamlss} fitting. Typically
     only used when calling \code{gamlss} function with the option \code{control}. 

}
\usage{
gamlss.control(c.crit = 0.001, n.cyc = 20, mu.step = 1, sigma.step = 1, nu.step = 1, 
               tau.step = 1, gd.tol = 5, iter = 0, trace = TRUE, autostep = TRUE, 
               save = TRUE, ...)

}

\arguments{
  \item{c.crit}{the convergence criterion for the algorithm}
  \item{n.cyc}{the number of cycles of the algorithm}
  \item{mu.step}{the step length for the parameter \code{mu}}
  \item{sigma.step}{the step length for the parameter \code{sigma} }
  \item{nu.step}{the step length for the parameter \code{nu} }
  \item{tau.step}{the step length for the parameter \code{tau} }
  \item{gd.tol}{global deviance tolerance level}
  \item{iter}{starting value for the number of iterations, typically set to 0 unless the function \code{\link{refit}} is used} 
  \item{trace}{whether to print at each iteration (TRUE) or not (FALSE) }
  \item{autostep}{whether the steps should be halved automatically if the new global deviance is greater that the old one, 
                   the default is \code{autostep=TRUE}  }
  \item{save}{\code{save=TRUE}, (the default), saves all the information on exit. 
              \code{save=FALSE} saves only limited information as the global deviance and AIC. 
              For example fitted values, design matrices and additive terms are not saved.
            The latest is useful when gamlss() is called several times within a procedure. }                  
  \item{\dots}{for extra arguments}
}
\details{
The step length for each of the parameters \code{mu}, \code{sigma}, \code{nu} or \code{tau} is very useful to aid convergence 
if the parameter has a fully parametric model.
However using a step length  is not theoretically justified if the model for the parameter includes one or more smoothing terms, 
(even thought it may give a very approximate result). 

The \code{c.crit} can be increased to speed up the convergence especially for a large set of data which takes longer to fit.   
When `trace' is TRUE, calls to the function \code{cat} produce the output for each outer iteration.
}
\value{
   A list with the arguments as components.
}
\references{ 
Rigby, R. A. and  Stasinopoulos D. M. (2005). Generalized additive models for location, scale and shape,(with discussion), 
\emph{Appl. Statist.}, \bold{54}, part 3, pp 507-554.

Stasinopoulos D. M., Rigby R.A. and Akantziliotou C. (2006) Instructions on how to use the GAMLSS package in R.
Accompanying documentation in the current GAMLSS  help files, (see also  \url{http://www.gamlss.org/}). 

Stasinopoulos D. M. Rigby R.A. (2007) Generalized additive models for location scale and shape (GAMLSS) in R.
\emph{Journal of Statistical Software}, Vol. \bold{23}, Issue 7, Dec 2007, \url{http://www.jstatsoft.org/v23/i07}.
}

\author{ Mikis Stasinopoulos \email{d.stasinopoulos@londonmet.ac.uk}, Bob Rigby \email{r.rigby@londonmet.ac.uk} 
}

\seealso{ \code{\link{gamlss}} }
\examples{
data(aids)
h<-gamlss(y~poly(x,3)+qrt, family=PO, data=aids) # 
con<-gamlss.control(mu.step=0.1)
h<-gamlss(y~poly(x,3)+qrt, family=PO, data=aids, control=con) # 
rm(h,con)
}
\keyword{regression}% 
