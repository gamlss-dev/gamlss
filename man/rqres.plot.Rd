\name{rqres.plot}
\alias{rqres.plot}
\alias{get.rqres}

\title{Creating and Plotting Randomized Quantile Residuals}
\description{
 This function plots worm plots, van Buuren and Fredriks M. (2001),  or QQ-plots of the normalized randomized quantile residuals (Dunn and Smyth, 1996) for a model using a discrete GAMLSS family distribution.    
}
\usage{
rqres.plot(obj = NULL, howmany = 6, plot.type = c("few", "all"), 
           type = c("wp", "QQ"),  xlim = NULL, ylim = NULL, ...)
get.rqres(obj = NULL, howmany = 10, order = FALSE)           
           
}

\arguments{
  \item{obj}{a fitted GAMLSS model object from a "discrete" type of family }
  \item{howmany}{The number randomise  quantile residuals required i.e. \code{howmany=6}}
   \item{plot.type}{whether to plot few of the randomised quantile residual realisations, \code{"few"} in a separate plots (there must be less than 8)  or all   \code{"all"} in one plot (with their median)}
     \item{type}{whether to plot worm plots \code{"wp"}or QQ plots \code{"QQ"} with default worm plots}
     \item{xlim}{setting manually the \code{xlim} of the graph}
     \item{ylim}{setting manually the \code{ylim} of the graph}
     \item{order}{whether to order the ealization of randomised quantile residuals}
   \item{\dots}{for extra arguments to be passed to \code{wp()}}
}
\details{For discrete family distributions, the \code{\link{gamlss}()} function saves on exit one realization of randomized quantile residuals which 
can be plotted using the generic function \code{plot} which calls the \code{plot.gamlss}. Looking at only one realization can be misleading, so the 
current function creates QQ-plots for several 
realizations. The function allows up to 10 QQ-plots to be plotted. Occasionally  one wishes to create a lot of realizations 
and then take a median of them (separately for each ordered value) to create a single median realization. The option \code{all} in combinations 
with the option \code{howmany} creates a 
QQ-plot of the medians of the normalized randomized quantile residuals. These 'median' randomized quantile residuals can be saved using the option
(\code{save=TRUE}).  
 
}
\value{
  
  If \code{save} it is TRUE then the vector of the median residuals is saved. 
}
\references{ Dunn, P. K. and Smyth, G. K. (1996) Randomised quantile residuals,
    \emph{J. Comput. Graph. Statist.}, \bold{5}, 236--244
    
Rigby, R. A. and  Stasinopoulos D. M. (2005). Generalized additive models for location, scale and shape,(with discussion), 
\emph{Appl. Statist.}, \bold{54}, part 3, pp 507-554.


Rigby, R. A., Stasinopoulos, D. M.,  Heller, G. Z.,  and De Bastiani, F. (2019)
	\emph{Distributions for modeling location, scale, and shape: Using GAMLSS in R}, Chapman and Hall/CRC. An older version can be found in \url{https://www.gamlss.com/}.

Stasinopoulos D. M. Rigby R.A. (2007) Generalized additive models for location scale and shape (GAMLSS) in R.
\emph{Journal of Statistical Software}, Vol. \bold{23}, Issue 7, Dec 2007, \url{https://www.jstatsoft.org/v23/i07/}.

Stasinopoulos D. M., Rigby R.A., Heller G., Voudouris V., and De Bastiani F., (2017)
\emph{Flexible Regression and Smoothing: Using GAMLSS in R},  Chapman and Hall/CRC.  

(see also \url{https://www.gamlss.com/}).

van Buuren and Fredriks M. (2001) Worm plot: simple diagnostic device for modelling growth reference curves. 
            \emph{Statistics in Medicine}, \bold{20}, 1259--1277
}

\author{Mikis Stasinopoulos \email{d.stasinopoulos@londonmet.ac.uk}}


\seealso{ \code{\link{plot.gamlss}}, \code{\link{gamlss}}  }

\examples{
data(aids) # fitting a model from a discrete distribution 
h<-gamlss(y~pb(x)+qrt, family=NBI, data=aids) # 
plot(h)
# plot qq- plots from 6 realization of the randomized quantile residuals
rqres.plot(h) 
# a worm-plot of the medians from 10 realizations
rqres.plot(h,howmany=40,plot="all") # 
}
\keyword{regression}% 
