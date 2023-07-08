\name{loglogSurv}
\alias{loglogSurv}
\alias{loglogSurv1}
\alias{loglogSurv2}
\alias{loglogSurv3}
\alias{logSurv}
\alias{logSurv0}
\alias{ECDF}
\alias{loglogplot}
\alias{loglogplot0}
\title{
Survival function plots for checking the tail behaviour of the data
}
\description{
The log-log Survival functions are design for checking the tails of a single response variable (no explanatory should be involved). There are three different function:

a) the function \code{loglogSurv1()} which plot the right tails of the empirical log-log Survival function against \code{loglog(y)}, where y is the variable of interest. The coefficient of a linear fit to the plot can be used as an estimated for Type I tails  (see Chapter 17 in Rigby \emph{et al.} (2019) for definition of the different types of tails.) 

b)  the function \code{loglogSurv2()} which plot the  right tails of the empirical log-log Survival function against \code{log(y)}. The coefficient of a linear fit to the plot can be used as an  estimated for Type II tails. 

c)  the function \code{loglogSurv3()} which plot the (left or right) tails of the empirical log-log Survival function against
\code{y}. The coefficient of a linear fit to the plot can be used as an estimated for Type III tails. 

The function \code{loglogSurv()} combines all the above functions. 

The function \code{logSurv()} is design for exploring the heavy tails of a single response variable. It plots the empirical log-survival function of the right tail of the distribution or the empirical log-cdf function  of the  left tail against \code{log(y)} for a specified probability of the tail. Then fits a linear, a quadratic and an exponential curve to the points of the plot. For distributions defined on the positive real line a good linear fit would indicate a Pareto type tail, a good quadratic fit a log-normal type tail and good exponential fit a Weibull type tail. Note that this function is only appropriate to investigate rather heavy tails and it is not very good to discriminate between different type of tails, as the \code{loglogSurv()}. The function \code{logSurv0()} plots but do not fit the curves.

The function \code{loglogplot()} plot  the empirical log-survival function of all data against  \code{log(y)}. 
The function \code{ECDF()} calculates the empirical commutative distribution function. It is similar to \code{ecdf()} but divides by \code{n+1} rather \code{n}, the number of conservations.   
}

\usage{
loglogSurv(y, prob = 0.9, print = TRUE, title = NULL, lcol = gray(0.1), 
           ltype = 1, plot = TRUE, ...)

loglogSurv1(y, prob = 0.9, print = TRUE, title = NULL, lcol = gray(0.1), 
           ltype = 1, ...)

loglogSurv2(y, prob = 0.9, print = TRUE, title = NULL, lcol = gray(0.1), 
           ltype = 1, ...)
           
loglogSurv3(y, prob = 0.9, print = TRUE, title = NULL, lcol = gray(0.1), 
           ltype = 1, ...)
          
logSurv(y, prob = 0.9, tail = c("right", "left"), plot = TRUE, 
       lines = TRUE, print = TRUE, title = NULL, lcol = c(gray(0.1), 
       gray(0.2), gray(0.3)), ltype = c(1, 2, 3), ...)  

logSurv0(y, prob = 0.9, tail = c("right", "left"), plot = TRUE, 
           title = NULL, ...)
           
         
ECDF(y, weights=NULL)

loglogplot(y, nplus1 = TRUE, ...)

}

\arguments{
  \item{y}{a vector, the variable of interest}
  \item{prob}{what probability. The defaul is  0.90 which means 10\% for "right" tail 90\% for "left" tail }
  \item{tail}{which tall needs checking the right (default) of the left}
  \item{plot}{whether to plot with default equal \code{TRUE} }
  \item{print}{whether to print the coefficients with default equal \code{TRUE}}
  \item{title}{if a different title rather the default is needed}
  \item{lcol}{The line colour in the plot}
  \item{lines}{whether to plot the fitted lines}  
  \item{ltype}{The line type in the plot}
  \item{nplus1}{whether to divide by n+1 or n when calculating the ecdf}
  \item{weights}{prior weights for \code{ECDF()}}  
  \item{...}{for extra argument in the plot command}
}
\details{


The functions \code{loglogSurv1()}, \code{loglogSurv3()} and \code{loglogSurv3()} take the upper  part of an ordered variable, create its empirical survival function, and plot the log-log of this functions against \code{log(log(y))}, \code{log(y)} and \code{y}, respectively. Then they fit a line to the plot. The coefficients of the line can be interpreted as parameters determined the behaviour of the tail.
The function \code{loglogSurv()} fits all three models and displays  the best. 


The function \code{logSurv()} takes the upper (or lower) part of an ordered variable  and plots the log empirical survival function against log(y). Also display three curves i) linear ii) quadratic and iii) exponential to determine what kind of tail relationship exist. Plotting the log  empirical survival function against log(y) often call in the literature the "log-log plot". 

The function \code{loglogplot()} plots the whole log  empirical survival function against log(y) (not just the tail).  The function \code{ECDF()}  calculate the step function of the empirical cumulative distribution function.


More details can be found in Chapter 17 of  "Rigby \emph{et al.} (2019)  book  an old version on which can be found in   \url{https://www.gamlss.com/})
}
\value{
The functions create plots. 
}
\references{
Rigby, R. A. and  Stasinopoulos D. M. (2005). Generalized additive models for location, scale and shape,(with discussion), 
\emph{Appl. Statist.}, \bold{54}, part 3, pp 507-554.


Rigby R.A., Stasinopoulos D. M., Heller G., and De Bastiani F., (2019) \emph{Distributions for Modelling Location, Scale and Shape: Using GAMLSS in R}, Chapman and Hall/CRC. (In press)

Rigby, R. A., Stasinopoulos, D. M.,  Heller, G. Z.,  and De Bastiani, F. (2019)
	\emph{Distributions for modeling location, scale, and shape: Using GAMLSS in R}, Chapman and Hall/CRC. An older version can be found in \url{https://www.gamlss.com/}.

Stasinopoulos D. M. Rigby R.A. (2007) Generalized additive models for location scale and shape (GAMLSS) in R.
\emph{Journal of Statistical Software}, Vol. \bold{23}, Issue 7, Dec 2007, \url{https://www.jstatsoft.org/v23/i07/}.

Stasinopoulos D. M., Rigby R.A., Heller G., Voudouris V., and De Bastiani F., (2017)
\emph{Flexible Regression and Smoothing: Using GAMLSS in R},  Chapman and Hall/CRC.  

(see also \url{https://www.gamlss.com/}).
}
\author{
Bob Rigby, Mikis Stasinopoulos and Vlassios Voudouris 
}


%\seealso{}
\examples{
data(film90)
y <- film90$lborev1
op<-par(mfrow=c(3,1))
loglogSurv1(y)
loglogSurv2(y)
loglogSurv3(y)
par(op)
loglogSurv(y)

logSurv(y)

loglogplot(y)

plot(ECDF(y), main="ECDF")
}

\keyword{distribution}

