\name{centiles}
\alias{centiles}
\alias{centiles.fan}

\title{Plots the centile curves for a GAMLSS object }
\description{This function \code{centiles()} plots centiles curves for distributions belonging to the GAMLSS family of distributions.  
The function also tabulates the sample percentages below each centile curve (for comparison with the model percentages given by the argument \code{cent}.)
The function \code{centiles.fan()} plots a fan-chart of the centile curves.
A restriction of the functions is that it applies to models with one explanatory variable only.    
}
\usage{
centiles(obj, xvar = NULL, cent = c(0.4, 2, 10, 25, 50, 75, 90, 98, 99.6), 
         legend = TRUE, ylab = "y", xlab = "x", main = NULL, 
         main.gsub = "@", xleg = min(xvar), yleg = max(obj$y), 
         xlim = range(xvar), ylim = range(obj$y), save = FALSE, 
         plot = TRUE, points = TRUE, pch = "+", col = "blue", 
         col.centiles = 1:length(cent) + 2, lty.centiles = 1, lwd.centiles = 1, ...)
centiles.fan(obj, xvar = NULL, cent = c(0.4, 2, 10, 25, 50, 75, 90, 98, 99.6), 
         ylab = "y", xlab = "x", main = NULL, main.gsub = "@", 
         xleg = min(xvar), yleg = max(obj$y), xlim = range(xvar), 
         ylim = range(obj$y), points = FALSE,  median = TRUE,  pch = "+", 
         col = "blue", 
         colors = c("cm", "gray", "rainbow", "heat", "terrain", "topo"), ...)
}

\arguments{
  \item{obj}{a fitted gamlss object from fitting a gamlss distribution}
  \item{xvar}{the unique explanatory variable}
  \item{cent}{a vector with elements the \% centile values for which the centile curves have to be evaluated}
  \item{legend}{whether a legend is required in the plot or not, the default is \code{legent=TRUE}  }
  \item{ylab}{the y-variable label}
  \item{xlab}{the x-variable label}
  \item{main}{the main title here as character. 
              If NULL the default title "centile curves using  NO" (or the relevant distributions name) is shown}
  \item{main.gsub}{if the \code{main.gsub} (with default "@") appears in the \code{main} 
                 title then it is substituted with the default title.}
  \item{xleg}{position of the legend in the x-axis}
  \item{yleg}{position of the legend in the y-axis}
  \item{xlim}{the limits of the x-axis}
  \item{ylim}{the limits of the y-axis}
  \item{save}{whether to save the sample percentages or not with default equal to \code{FALSE}. 
            In this case the sample percentages are printed but are not saved }
  \item{plot}{whether to plot the centiles. This option is useful for \code{centile.split} } 
  \item{pch}{the character to be used as the default in plotting points see \code{par}}
  \item{col}{plotting colour see \code{par}}
  \item{col.centiles}{Plotting colours  for the centile curves}
  \item{lty.centiles}{line type for the centile curves}
  \item{lwd.centiles}{The line width for the centile curves}
  \item{colors}{the different colour schemes to be used for the fan-chart. The following are available 
  \code{c("cm","gray", "rainbow", "heat", "terrain", "topo")},  }
  \item{points}{whether the data points should be plotted, default is \code{TRUE} for \code{centiles()} and
                \code{FALSE} for \code{centiles.fan()}}
  \item{median}{whether the median should be plotted (only in \code{centiles.fan()})}
  \item{\dots}{for extra arguments}          
}
\details{
 Centiles are calculated using the fitted values in \code{obj} and \code{xvar} must
correspond exactly to the predictor in \code{obj} to plot correctly.

\code{col.centiles}, \code{lty.centiles} and \code{lwd.centiles} may be vector arguments
and are recycled to the length \code{cent} if necessary.
}
\value{
A centile plot is produced and the sample centiles below each centile curve are printed (or saved)
}
\references{Rigby, R. A. and  Stasinopoulos D. M. (2005). Generalized additive models for location, scale and shape,(with discussion), 
\emph{Appl. Statist.}, \bold{54}, part 3, pp 507-554.

Stasinopoulos D. M., Rigby R.A. and Akantziliotou C. (2006) Instructions on how to use the GAMLSS package in R.
Accompanying documentation in the current GAMLSS  help files, (see also  \url{http://www.gamlss.org/}). 

Stasinopoulos D. M. Rigby R.A. (2007) Generalized additive models for location scale and shape (GAMLSS) in R.
\emph{Journal of Statistical Software}, Vol. \bold{23}, Issue 7, Dec 2007, \url{http://www.jstatsoft.org/v23/i07}.
}

\author{Mikis Stasinopoulos \email{d.stasinopoulos@londonmet.ac.uk}, Bob Rigby \email{r.rigby@londonmet.ac.uk} with contribution from Steve Ellison}


\section{Warning }{This function is appropriate only when one continuous explanatory variable is fitted in the model}  
\seealso{  \code{\link{gamlss}}, \code{\link{centiles.split}} , \code{\link{centiles.com}}}
\examples{
data(abdom)
h<-gamlss(y~pb(x), sigma.formula=~pb(x), family=BCT, data=abdom) 
# default plot
centiles(h,xvar=abdom$x)
# control of colours and lines
centiles(h, xvar=abdom$x,  col.cent=c(2,3,4,5,1,5,4,3,2,1), 
              lwd.cent=c(1,1,1,1,2,1,1,1,1))
#Control line types
centiles(h, xvar=abdom$x,  col.cent=1, cent=c(.5,2.5,50,97.5,99.5), 
              lty.centiles=c(3,2,1,2,3),lwd.cent=c(1,1,2,1,1))
# control of the main title
centiles(h, xvar=abdom$x,  main="Abdominal data \n @")
# the fan-chart
centiles.fan(h,xvar=abdom$x, colors="rainbow")
rm(h)
}
\keyword{regression}% 
