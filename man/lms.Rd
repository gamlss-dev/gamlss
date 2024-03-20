\name{lms}
\alias{lms}

\title{A function to fit LMS curves for centile estimation}
\description{
This function is design to help the user to easily construct growth curve centile estimation.
It is  applicable when only "one" explanatory variable is available (usually age).
}
\usage{
lms(y, x, families = LMS, data = NULL, k = 2, 
        cent = c(0.4, 2, 10, 25, 50, 75, 90, 98, 99.6),
        calibration = TRUE, trans.x = FALSE, 
        fix.power = NULL, lim.trans = c(0, 1.5), 
        prof = FALSE, step = 0.1, legend = FALSE, 
        mu.df = NULL, sigma.df = NULL, nu.df = NULL, 
        tau.df = NULL, c.crit = 0.01, 
        method.pb = c("ML", "GAIC"), ...)        
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{The response variable}
  \item{x}{The unique explanatory variable}
  \item{families}{a list of \code{gamlss.families} with default \code{LMS=c("BCCGo",  "BCPEo", "BCTo")}}
  \item{data}{the data frame}
  \item{k}{the penalty to be used in the GAIC}
 \item{cent}{a vector with elements the \% centile values for which the centile curves have to be evaluated} 
 \item{calibration}{whether calibration is required with default \code{TRUE}}
 \item{trans.x}{whether to check for transformation in x with default \code{FALSE}}
  \item{fix.power}{if set it fix the power of the transformation for x}
 \item{lim.trans}{the limits for the search of the power parameter for x}
  \item{prof}{whether to use the profile GAIC of the power transformation}
  \item{step}{if \code{prof=TRUE} is used this determine the step for the profile GAIC }
 \item{legend}{whether a legend is required in the plot with default \code{FALSE}}
 \item{mu.df}{\code{mu} effective degrees of freedom if required otherwise are estimated}
  \item{sigma.df}{\code{sigma} effective degrees of freedom if required otherwise are estimated}
  \item{nu.df}{\code{nu} effective degrees of freedom if required otherwise are estimated}
  \item{tau.df}{\code{tau} effective degrees of freedom if required otherwise are estimated}
   \item{c.crit}{the convergence criterion to be pass to \code{gamlss()} }
  \item{method.pb}{the method used in the \code{pb()} for estimating the smoothing 
  parameters. The default is local maximum likelihood \code{"ML"}. \code{"GAIC"} is also permitted where \code{k} is taken from the \code{k} argument of the function.}
  \item{\dots}{extra argument which can be passed to \code{gamlss()}}
}
\details{
This function should be used if the construction of the centile curves involves only one explanatory variable. 

The model assumes that the response variable has a flexible distribution i.e. \eqn{y ~ D(\mu, \sigma, \nu, \tau)} where the parameters of the distribution are smooth functions of the explanatory variable i.e. \eqn{g(\mu)= s(x)}, where \eqn{g()} is a link function and \eqn{s()} is a smooth function. Occasionally a power transformation in the x-axis helps the construction of the centile curves. That is, in this case the parameters are modelled by \eqn{x^p} rather than just x, i.e.\eqn{g(\mu)= s(x^p)}. The function \code{lms()} uses  P-splines (\code{pb()}) as a smoother. 

If a transformation is needed for \code{x} the function \code{lms()}  starts by finding an optimum value for \code{p}  using the simple model  \eqn{NO(\mu=s(x^p))}. (Note that this value of \code{p} is not the optimum for the final chosen model but it works well in practice.) 

After fitting a Normal error model for staring values the function proceeds by fitting several "appropriate" distributions for the response variable.
The set of \code{gamlss.family} distributions to fit is specified by the argument \code{families}. 
The default \code{families} arguments is \code{LMS=c("BCCGo",  "BCPEo", "BCTo")} that is the LMS class of distributions, Cole and Green (1992).
Note that this class is only appropriate when y is positive (with no zeros). If the response variable contains negative values and zeros then use the argument \code{families=theSHASH} where \code{theSHASH <-  c("NO", "SHASHo")} or add any other list of distributions which you may think is appropriate.
Justification of using the specific centile (0.38  2.27  9.1211220 25.25, 50, 74.75, 90.88, 97.72, 99.62) is given in Cole (1994).
}
\value{It returns a \code{gamlss} fitted object}

\references{
Cole, T. J. (1994) Do growth chart centiles need a face lift?
 \emph{BMJ}, 308--641.

Cole,  T. J. and Green,  P. J. (1992) Smoothing reference centile curves: the LMS method and penalized likelihood, \emph{Statist. Med.} \bold{11}, 1305--1319

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
\author{
Mikis Stasinopoulos, Bob Rigby and Vlasios Voudouris}

\note{The function is fitting several models and for large data can be slow}

\seealso{\code{\link{gamlss}}, \code{\link{centiles}}, \code{\link{calibration}}}
\examples{
\dontrun{
data(abdom)
m1 <- lms(y,x , data=abdom, n.cyc=30)
m2 <- lms(y,x ,data=abdom, method.pb="GAIC", k=log(610))
# this example takes time
data(db)
m1 <- lms(y=head, x=age, data=db, trans.x=TRUE)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{regression}
