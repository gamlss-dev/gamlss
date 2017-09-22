\name{term.plot}
\alias{term.plot}

\title{Plot regression terms for a specified parameter of a fitted GAMLSS object}
\description{
  Plots regression terms against their predictors, optionally with
     standard errors and partial residuals added. It is based on the R function \code{termplot} but is suitably changed to apply to GAMLSS objects.

}
\usage{
term.plot(object, what = c("mu", "sigma", "nu", "tau"), 
         parameter= NULL, data = NULL, 
         envir = environment(formula(object)), partial.resid = FALSE, 
         rug = FALSE, terms = NULL, se = TRUE, ylim = c("common", "free"), 
         scheme = c("shaded", "lines"), xlabs = NULL, ylabs = NULL, 
         main = NULL, pages = 0, col.term = "darkred", 
         col.se = "orange", col.shaded = "gray", col.res = "lightblue", 
         col.rug = "gray", lwd.term = 1.5, lty.se = 2, lwd.se = 1, 
         cex.res = 1, pch.res = par("pch"), 
         ask = interactive() && nb.fig < n.tms && .Device != "postscript", 
         use.factor.levels = TRUE, surface.gam = FALSE,     
         polys = NULL, polys.scheme = "topo",...)

}

\arguments{
  \item{object}{a fitted GAMLSS object}
  \item{what}{the required parameter of the GAMLSS distribution i.e. "mu"}
  \item{parameter}{equivalent to \code{what}}
  \item{data}{data frame in which variables in \code{object} can be found}
  \item{envir}{environment in which variables in \code{object} can be found }
  \item{partial.resid}{logical; should partial residuals be plotted or not}
  \item{rug}{ add rug plots (jitter 1-d histograms) to the axes?}
  \item{terms}{which terms to be plotted (default 'NULL' means all terms) }
  \item{se}{plot point-wise standard errors?}
  \item{ylim}{there are two options here a) "common" and b) "free". 
       The "common" option plots all figures with the same \code{ylim} 
       range and therefore allows the viewer to check the relative 
       contribution of each terms compate to the rest.  
       In the`free' option the limits are computed for each plot seperatly.}
  \item{scheme}{whether the se's should appear shaded or as lines}
  \item{xlabs}{vector of labels for the x axes }
  \item{ylabs}{vector of labels for the y axes }
  \item{main}{logical, or vector of main titles;  if 'TRUE', the model's
          call is taken as main title, 'NULL' or 'FALSE' mean no
          titles.}
   \item{pages}{in how many pages the plot should appear. The default 
            is 0 which allows differnt page for each plot}
  \item{col.term}{the colour of the term line}
  \item{col.se}{the colour of the se's lines}
  \item{col.shaded}{the colour of the shaded area}
  \item{col.res}{the colour of the partial residuals}
  \item{col.rug}{the colour of the rug}
  \item{lwd.term}{line width of the fitted terms}
  \item{lty.se}{line ype for standard errors}
  \item{lwd.se}{line width for the stadard errors}
  \item{cex.res}{ plotting character expansion for the parsial residuals}
  \item{pch.res}{characters for points in the parsial residuals}
  \item{ask}{logical; if 'TRUE', the user is asked before each plot, see
          'par(ask=.)'.}        
  \item{use.factor.levels}{Should x-axis ticks use factor levels or numbers for
          factor terms? }
  \item{surface.gam}{whether to use surface plot if a \code{ga()} term is fitted}
  \item{polys}{The polygone nformation filr for MRF models}
   \item{polys.scheme}{Color scheme for polygones for RMF models}
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

Stasinopoulos D. M. Rigby R.A. (2007) Generalized additive models for location scale and shape (GAMLSS) in R.
\emph{Journal of Statistical Software}, Vol. \bold{23}, Issue 7, Dec 2007, \url{http://www.jstatsoft.org/v23/i07}.

Stasinopoulos D. M., Rigby R.A., Heller G., Voudouris V., and De Bastiani F., (2017) \emph{Flexible Regression and Smoothing: Using GAMLSS in R}, Chapman and Hall/CRC. 

(see also \url{http://www.gamlss.org/}).
}
\author{Mikis Stasinopoulos based on the existing termplot() function}

\seealso{ \code{\link{termplot}}}
\examples{
data(aids)
a<-gamlss(y~pb(x)+qrt,data=aids,family=NBI)
term.plot(a, pages=1)
rm(a)
}
\keyword{regression}% 
