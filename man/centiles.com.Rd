\name{centiles.com}
\alias{centiles.com}

\title{
Comparing centiles from different GAMLSS models}
\description{
This function compares centiles curves for more than one GAMLSS objects.It is based on the \code{centiles} function.
 The function also tabulates the sample percentages below each centile curve (for comparison with the model percentages 
given by the argument \code{cent}.) A restriction of the function is that it applies to models with one 
explanatory variable only
}

\usage{
centiles.com(obj, ..., xvar, cent = c(0.4, 10, 50, 90, 99.6), 
             legend = TRUE, ylab = "y", xlab = "x", xleg = min(xvar), 
             yleg = max(obj$y), xlim = range(xvar), ylim = NULL, 
             no.data = FALSE, color = TRUE, main = NULL, plot = TRUE)
}

\arguments{
  \item{obj}{a fitted gamlss object from fitting a gamlss continuous distribution}
  \item{\dots}{optionally more fitted GAMLSS model objects}
  \item{xvar}{the unique explanatory variable }
  \item{cent}{a vector with elements the \% centile values for which the centile curves have to be evaluated}
  \item{legend}{whether a legend is required in the plot or not, the default is \code{legent=TRUE}}
  \item{ylab}{the y-variable label}
  \item{xlab}{the x-variable label}
  \item{xleg}{position of the legend in the x-axis}
  \item{yleg}{position of the legend in the y-axis}
  \item{xlim}{the limits of the x-axis}
  \item{ylim}{the limits of the y-axis}
  \item{no.data}{whether the data should plotted, default \code{no.data=FALSE} or not \code{no.data=TRUE} }
  \item{color}{whether the fitted centiles are shown in colour, \code{color=TRUE} (the default) or not          
               \code{color=FALSE}}
  \item{main}{the main title}
  \item{plot}{whether to plot the centiles}
  }


\value{Centile plots are produced for the different fitted models and the sample centiles below each centile curve are printed}
\references{
Rigby, R. A. and  Stasinopoulos D. M.(2005). Generalized additive models for location, scale and shape, (with discussion),\emph{Appl. Statist.}, 
\bold{54}, part 3, pp 507-554.

Rigby, R. A., Stasinopoulos, D. M.,  Heller, G. Z.,  and De Bastiani, F. (2019)
	\emph{Distributions for modeling location, scale, and shape: Using GAMLSS in R}, Chapman and Hall/CRC. An older version can be found in \url{https://www.gamlss.com/}.

Stasinopoulos D. M. Rigby R.A. (2007) Generalized additive models for location scale and shape (GAMLSS) in R.
\emph{Journal of Statistical Software}, Vol. \bold{23}, Issue 7, Dec 2007, \url{https://www.jstatsoft.org/v23/i07/}.

Stasinopoulos D. M., Rigby R.A., Heller G., Voudouris V., and De Bastiani F., (2017)
\emph{Flexible Regression and Smoothing: Using GAMLSS in R},  Chapman and Hall/CRC.  

(see also \url{https://www.gamlss.com/}).

}
\author{
Mikis Stasinopoulos and Bob Rigby}

\section{Warning }{This function is appropriate only when one continuous explanatory variable is fitted in the model} 

\seealso{ \code{\link{gamlss}}, \code{\link{centiles}} , \code{\link{centiles.split}}}
\examples{
data(abdom)
h1<-gamlss(y~cs(x,df=3), sigma.formula=~cs(x,1),family=BCT, data=abdom)
h2<-gamlss(y~pb(x), sigma.formula=~pb(x), family=BCT, data=abdom )
centiles.com(h1,h2,xvar=abdom$x)
rm(h1,h2)
}
\keyword{regression}
