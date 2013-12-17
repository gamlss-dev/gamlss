\name{term.plot}
\alias{term.plot}

\title{Plot regression terms for a specified parameter of a GAMLSS object}
\description{
  Plots regression terms against their predictors, optionally with
     standard errors and partial residuals added. It is almost identical to the R function \code{termplot} 
     suitably changed to apply to GAMLSS objects.

}
\usage{
term.plot(object, what = c("mu", "sigma", "nu", "tau"), data = NULL, 
          envir = environment(formula(object)), 
          partial.resid = FALSE, rug = FALSE, 
          terms = NULL, se = TRUE, xlabs = NULL, ylabs = NULL, 
          main = NULL, col.term = 2, lwd.term = 1.5, 
          col.se = "orange", lty.se = 2, lwd.se = 1, 
          col.res = "gray", cex = 1, pch = par("pch"), 
          col.smth = "darkred", lty.smth = 2, 
          span.smth = 2/3, 
          ask = interactive() && nb.fig < n.tms && .Device != "postscript", 
          use.factor.levels = TRUE, smooth = NULL,  ylim  = "common", ...)
}

\arguments{
  \item{object}{a GAMLSS object}
  \item{what}{the required parameter of the GAMLSS distribution}
  \item{data}{data frame in which variables in \code{object} can be found}
  \item{envir}{environment in which variables in \code{object} can be found }
  \item{partial.resid}{logical; should partial residuals be plotted?}
  \item{rug}{ add rug plots (jitter 1-d histograms) to the axes?}
  \item{terms}{which terms to plot (default 'NULL' means all terms) }
  \item{se}{plot point-wise standard errors?}
  \item{xlabs}{vector of labels for the x axes }
  \item{ylabs}{vector of labels for the y axes }
  \item{main}{logical, or vector of main titles;  if 'TRUE', the model's
          call is taken as main title, 'NULL' or 'FALSE' mean no
          titles.}
  \item{col.term, lwd.term}{colour and width for the "term curve", see 'lines'.}
  \item{col.se, lty.se, lwd.se}{ colour, line type and line width for the
          "twice-standard-error curve" when 'se = TRUE'.}
  \item{col.res, cex, pch}{colour, plotting character expansion and type for
          partial residuals, when 'partial.resid = TRUE', see 'points'.}
  \item{lty.smth,col.smth, span.smth}{Passed to 'smooth'}
  \item{ask}{logical; if 'TRUE', the user is asked before each plot, see
          'par(ask=.)'.}
  \item{use.factor.levels}{Should x-axis ticks use factor levels or numbers for
          factor terms? }
  \item{smooth}{'NULL' or a function with the same arguments as
          'panel.smooth' to draw a smooth through the partial residuals
          for non-factor terms}
  \item{ylim}{an optional range for the y axis, or "common" when a range sufficient for all the plot will be computed, or "free" when limits are computed for each plot.}
  \item{\dots}{other graphical parameters}
}
\details{
 The function uses the \code{lpred} function of GAMLSS.  
     The 'data' argument should rarely be needed, but in some cases
     'termplot' may be unable to reconstruct the original data frame.
     Using 'na.action=na.exclude' makes these problems less likely.
     Nothing sensible happens for interaction terms.
}
\value{
 a plot of fitted terms.
}
\references{ Rigby, R. A. and  Stasinopoulos D. M. (2005). Generalized additive models for location, scale and shape,(with discussion), 
\emph{Appl. Statist.}, \bold{54}, part 3, pp 507-554.

Stasinopoulos D. M., Rigby R.A. and Akantziliotou C. (2006) Instructions on how to use the GAMLSS package in R.
Accompanying documentation in the current GAMLSS  help files, (see also  \url{http://www.gamlss.org/}).

Stasinopoulos D. M. Rigby R.A. (2007) Generalized additive models for location scale and shape (GAMLSS) in R.
\emph{Journal of Statistical Software}, Vol. \bold{23}, Issue 7, Dec 2007, \url{http://www.jstatsoft.org/v23/i07}.
}
\author{Mikis Stasinopoulos based on the existing termplot() function}

\seealso{ \code{\link{termplot}}}
\examples{
data(aids)
a<-gamlss(y~pb(x)+qrt,data=aids,family=NBI)
term.plot(a, se=TRUE)
rm(a)
}
\keyword{regression}% 
