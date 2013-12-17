\name{pdf.plot}
\alias{pdf.plot}

\title{Plots Probability Distribution Functions for GAMLSS Family}

\description{A function to plot probability distribution functions (pdf) belonging to the gamlss family of distributions.
   This function allows either plotting of the fitted distributions for up to eight observations or plotting specified distributions belonging in the gamlss family}


\usage{
pdf.plot(obj = NULL, obs = c(1), family = NO(), mu = NULL, 
         sigma = NULL, nu = NULL, tau = NULL, min = NULL, 
         max = NULL, step = NULL, allinone = FALSE, 
         no.title = FALSE, ...)
}

\arguments{
  \item{obj}{An gamlss object e.g. \code{obj=model1} where model1 is a fitted gamlss object}
  \item{obs}{A number or vector of up to length eight indicating the case numbers of the observations for which fitted distributions are to be displayed, e.g. \code{obs=c(23,58)} will display the fitted distribution for the 23th and 58th observations}
  \item{family}{This must be a gamlss family i.e. \code{family=NO}}
  \item{mu}{The value(s) of the location parameter mu for which the distribution has to be evaluated e.g \code{mu=c(3,7)}}
  \item{sigma}{The value(s) the scale parameter sigma for which the distribution has to be evaluated e.g \code{sigma=c(3,7)}}
  \item{nu}{The value(s) the parameter nu for which the distribution has to be evaluated e.g. \code{nu=3}}
  \item{tau}{The value(s) the parameter tau for which the distribution has be evaluated e.g. \code{tau=5}}
  \item{min}{Minimum value of the random variable y e.g. \code{min=0}}
  \item{max}{Maximum value of y e.g. \code{max=10}}
  \item{step}{Steps for the evaluation of y e.g. \code{step=0.5}}
  \item{allinone}{This will go}
  \item{no.title}{Whether you need title in the plot, default is  \code{no.title=FALSE}}
  \item{...}{for extra arguments}
  }

\details{
   This function can be used to plot distributions of the GAMLSS family. 
   If the first argument \code{obj} is specified and it is a GAMLSS fitted object, then the fitted distribution of this model 
   at specified observation values (given by the second argument \code{obs}) is  plotted for a specified y-variable range (arguments 
  \code{min}, \code{max}, and \code{step}).
  
  If the first argument is not given then the \code{family} argument has to be specified and the pdf is plotted at specified values of the parameters 
  \code{mu}, \code{sigma}, \code{nu}, \code{tau}. Again the range of the y-variable has to be given.   
}

\value{plot(s) of the required pdf(s) are returned }

\references{
Rigby, R. A. and  Stasinopoulos D. M. (2005). Generalized additive models for location, scale and shape,(with discussion), 
\emph{Appl. Statist.}, \bold{54}, part 3, pp 507-554.

Stasinopoulos D. M., Rigby R.A. and Akantziliotou C. (2006) Instructions on how to use the GAMLSS package in R.
Accompanying documentation in the current GAMLSS  help files, (see also  \url{http://www.gamlss.org/}).

Stasinopoulos D. M. Rigby R.A. (2007) Generalized additive models for location scale and shape (GAMLSS) in R.
\emph{Journal of Statistical Software}, Vol. \bold{23}, Issue 7, Dec 2007, \url{http://www.jstatsoft.org/v23/i07}.
}
\author{Mikis Stasinopoulos \email{d.stasinopoulos@londonmet.ac.uk} and Calliope Akantziliotou}
\note{The range of the y values given by min, max and step are very important in the plot}

\section{Warning}{The range of some distributions depends on the fitted parameters} 
\seealso{ \code{\link[gamlss]{gamlss}} }
\examples{
pdf.plot(family=BCT, min=1, max=20, step=.05, mu=10, sigma=0.15, nu=-1, tau=c(4,10,20,40) )
# now using an gamlss object
# library(gamlss)
#data(abdom)
#h<-gamlss(y~pb(x), sigma.formula=~pb(x), family=BCT, data=abdom) # fits 
#pdf.plot(obj=h , obs=c(23,67), min=50, max=150, step=.5)
}
\keyword{regression}% 
