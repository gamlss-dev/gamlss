\name{fittedPlot}
\alias{fittedPlot}

\title{Plots The Fitted Values of a GAMLSS Model}
\description{
 This function, applicable only to a models with a single explanatory variable, plots the fitted values for all the parameters
  of a GAMLSS model against the (one) explanatory variable. It is also useful for comparing the fits for more than one model.  
}
\usage{
fittedPlot(object, ..., x = NULL, color = TRUE, line.type = FALSE, xlab = NULL)
}

\arguments{
  \item{object}{a fitted GAMLSS model object(with only one explanatory variable)}
  \item{\dots}{optionally more fitted GAMLSS model objects }
  \item{x}{The unique explanatory variable }
  \item{color}{whether the fitted lines plots are shown in colour, \code{color=TRUE} (the default) or not \code{color=FALSE}}
  \item{line.type}{whether the line type should be different or not. The default is  \code{color=FALSE}}
  \item{xlab}{the x-label}
}

\value{
 A plot of the fitted values against the explanatory variable 
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
\author{Mikis Stasinopoulos, Bob Rigby and Calliope Akantziliotou }

\seealso{  \code{\link{gamlss}}, \code{\link{centiles}}, \code{\link{centiles.split}} }
\examples{
data(abdom)
h1<-gamlss(y~pb(x), sigma.formula=~x, family=BCT, data=abdom)
h2<-gamlss(y~pb(x), sigma.formula=~pb(x), family=BCT, data=abdom)
fittedPlot(h1,h2,x=abdom$x)
rm(h1,h2)
}
\keyword{regression}% 
