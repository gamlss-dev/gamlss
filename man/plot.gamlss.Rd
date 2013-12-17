\name{plot.gamlss}
\alias{plot.gamlss}

\title{Plot Residual Diagnostics for an GAMLSS Object}
\description{
  This function provides four plots for checking the normalized (randomized for a discrete response distribution) quantile 
  residuals of a fitted GAMLSS object, referred to as residuals below : a plot of residuals against fitted values, a plot of the residuals against 
  an index or a specific explanatory variable,  a density plot of the residuals and a normal Q-Q plot of the residuals.
  If argument \code{ts=TRUE} then the first two plots are replaced by the autocorrelation function (ACF) and partial autocorrelation function (PACF)
  of the residuals 
}
\usage{
\method{plot}{gamlss}(x, xvar = NULL, parameters = NULL, ts = FALSE, 
           summaries = TRUE, ...)
}

\arguments{
  \item{x}{a GAMLSS fitted object }
  \item{xvar}{an explanatory variable to plot the residuals against}
  \item{parameters}{plotting parameters can be specified here}
  \item{ts}{set this to TRUE if ACF and PACF plots of the residuals are required }
  \item{summaries}{set this to FALSE if no summary statistics of the residuals are required}
  \item{...}{further arguments passed to or from other methods.}
}
\details{
This function provides four plots for checking the normalized (randomized) quantile residuals (called \code{residuals}) of a fitted GAMLSS object. 
Randomization is only performed for discrete response variables.  The four plots are 

\itemize{ 
 \item residuals against the fitted values (or ACF of the residuals if \code{ts=TRUE})
 \item residuals against an index or specified x-variable (or PACF of the residuals if \code{ts=TRUE})
 \item kernel density estimate of the residuals 
 \item QQ-normal plot of the residuals 
 }
For time series response variables option \code{ts=TRUE} can be used to plot the ACF and PACF functions of the residuals. 
}
\value{
   Returns four plots related to the residuals of the fitted GAMLSS model and prints summary statistics for the residuals if the \code{summary=T}
}
\references{ 

Rigby, R. A. and  Stasinopoulos D. M. (2005). Generalized additive models for location, scale and shape,(with discussion), 
\emph{Appl. Statist.}, \bold{54}, part 3, pp 507-554.

Stasinopoulos D. M., Rigby R.A. and Akantziliotou C. (2006) Instructions on how to use the GAMLSS package in R.
Accompanying documentation in the current GAMLSS  help files, (see also  \url{http://www.gamlss.org/}).
 
Stasinopoulos D. M. Rigby R.A. (2007) Generalized additive models for location scale and shape (GAMLSS) in R.
\emph{Journal of Statistical Software}, Vol. \bold{23}, Issue 7, Dec 2007, \url{http://www.jstatsoft.org/v23/i07}. 
 }
\author{ Mikis Stasinopoulos \email{d.stasinopoulos@londonmet.ac.uk}, Bob Rigby \email{r.rigby@londonmet.ac.uk} and Kalliope Akantziliotou }


\seealso{ \code{\link{gamlss}} }

\examples{

data(aids)
a<-gamlss(y~pb(x)+qrt,family=PO,data=aids)
plot(a)
rm(a)
}
\keyword{regression}% 
