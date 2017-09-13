\name{loglogSurv}
\alias{loglogSurv}
\alias{loglogSurv1}
\alias{loglogSurv2}
\alias{loglogSurv3}
\alias{logSurv}


\title{
Log-Log Survival function plots for checking the tail behaviour of the data
}
\description{
The log-log Survival functions are design for checking the tails of a single response variable (no explanatory should be involved). There are three different function:

a) the function \code{loglogSurv1()} which plot the (left or right) tails of the empirical log-log Survival function against
loglog(y), where y is the variable of interest. The coefficient of a linear fit to the plot can be used an estimated for  Type I tails. 

b)   the function \code{loglogSurv2()} which plot the (left or right) tails of the empirical log-log Survival function against
log(y). The coefficient of a linear fit to the plot can be used an  estimated for  Type II tails. 

c)  the function \code{loglogSurv3()} which plot the (left or right) tails of the empirical log-log Survival function against
y. The coefficient of a linear fit to the plot can be used an an estimated for  Type III tails. 

The function \code{loglogSurv()} combines all the above functions. 

The function \code{logSurv()} is also design for exploring the tails of a single response variable. It plots the empirical log-survival function against \code{log(y)} for specified percentage of the tail and fits a linear, quadratic and exponential curve to the points of the plot. For distributions defined on the positive real line a good linear fit would indicate a Pareto type tail, a good quadratic fit a log-normal type tail and good exponential fit a Weibull type tail. Note that this function is only appropriate to investigate rather heavy tails and it is not very good to discriminate between different type of tails, as the \code{loglogSurv()} .
}

\usage{
loglogSurv(y, percentage = 10, howmany = NULL, type = c("right", "left"), 
           plot = TRUE, print = TRUE, save = FALSE)
loglogSurv1(y, percentage = 10, howmany = NULL, type = c("right", "left"), 
           plot = TRUE, print = TRUE, save = FALSE)
loglogSurv2(y, percentage = 10, howmany = NULL, type = c("right", "left"), 
          plot = TRUE, print = TRUE, save = FALSE)
loglogSurv3(y, percentage = 10, howmany = NULL, type = c("right", "left"), 
          plot = TRUE, print = TRUE, save = FALSE)
          
logSurv(y, percentage = 10, howmany = NULL, type = c("right", "left"), 
          plot = TRUE, print = TRUE, save = FALSE)          
}

\arguments{
  \item{y}{a vector, the variable of interest}
  \item{percentage}{what percentage of the tail need to be modelled, default is 10\%}
  \item{howmany}{how many observations in the tail needed. This is an alternative to \code{percentage}. If it specified it take over from the \code{percentage} argument otherwise \code{percentage} is used.}
  \item{type}{which tall needs checking the right (default) of the left}
  \item{plot}{whether to plot with default equal \code{TRUE} }
  \item{print}{whether to print the coefficients with default equal \code{TRUE}}
  \item{save}{whether to save the fitted linear model with default equal \code{FALSE}}
}
\details{
The functions \code{loglogSurv1()}, \code{loglogSurv3()} and \code{loglogSurv3()} take the upper (or lower) part of an ordered variable create its empirical survival function and plot the log-log of this functions against \code{log(log(y))}, \code{log(y)} and \code{y} respectively. Then they fit a line to the plot. The coefficients of the line can be interpreted as parameters determined the behaviour of the tail. More details can be found in Chapter 6 of  "The Distribution Toolbox of GAMLSS" book  which can be found in   \url{http://www.gamlss.org/})
}
\value{
A plot
}
\references{
Rigby, R. A. and  Stasinopoulos D. M. (2005). Generalized additive models for location, scale and shape,(with discussion), 
\emph{Appl. Statist.}, \bold{54}, part 3, pp 507-554.

Stasinopoulos D. M. Rigby R.A. (2007) Generalized additive models for location scale and shape (GAMLSS) in R.
\emph{Journal of Statistical Software}, Vol. \bold{23}, Issue 7, Dec 2007, \url{http://www.jstatsoft.org/v23/i07}.

Stasinopoulos D. M., Rigby R.A., Heller G., Voudouris V., and De Bastiani F., (2017) \emph{Flexible Regression and Smoothing: Using GAMLSS in R}, Chapman and Hall/CRC. 

(see also \url{http://www.gamlss.org/}).
}
\author{
Bob Rigby, Mikis Stasinopoulos and Vlassios Voudouris 
}


%\seealso{}
\examples{
data(film90)
F90 <- film90$lborev1
op<-par(mfrow=c(3,1))
loglogSurv1(F90)
loglogSurv2(F90)
loglogSurv3(F90)
par(op)
loglogSurv(F90)

logSurv(F90)
}

\keyword{distribution}

