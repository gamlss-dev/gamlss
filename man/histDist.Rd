\name{histDist}
\alias{histDist}

\title{This function plots the histogram and a fitted (GAMLSS family) distribution to a variable}
\description{
  This function fits constants to the parameters of a GAMLSS family distribution and them plot the histogram and 
  the fitted distribution. 
}
\usage{
histDist(y, family = NO, freq = NULL, density = FALSE, 
        nbins = 10, xlim = NULL, ylim = NULL, main = NULL, 
        xlab = NULL, ylab = NULL, data = NULL, 
        col.hist = "gray", border.hist = "blue", 
        fg.hist = rainbow(12)[9], line.wd = 2, 
        line.ty = c(1, 2), line.col = c(2, 3), 
        col.main = "blue4", col.lab = "blue4", 
        col.axis = "blue", ...)
        }
\arguments{
  \item{y}{a vector for the response variable}
  \item{family}{a \code{gamlss.family} distribution}
  \item{freq}{the frequencies of the data in \code{y} if exist. \code{freq} is used as \code{weights} in the \code{gamlss} fit}
 \item{density}{default value is FALSE. Change to TRUE if you would like a non-parametric density plot together with the parametric 
                fitted distribution plot (for continuous variable only)}
  \item{nbins}{The suggested number of bins (argument passed to \code{truehist()} of package MASS). 
               Either a positive integer, or a character string naming a rule: "Scott" 
               or "Freedman-Diaconis" or "FD". (Case is ignored.) }
  \item{xlim}{the minimum and the maximum x-axis value (if the default values are out of range)}
  \item{ylim}{the minimum and the maximum y-axis value (if the default values are out of range)}
 \item{main}{the main title for the plot}
  \item{xlab}{the label in the x-axis}
  \item{ylab}{the label in the y-axis}
  \item{data}{the data.frame}
  \item{col.hist}{the colour of the histogram or barplot}
  \item{border.hist}{the colour of the border of the histogram or barplot}
  \item{fg.hist}{the colour of  axis in the histogram or barplot}
  \item{line.wd}{the line width of the fitted distribution}
  \item{line.ty}{the line type of the fitted distribution}
  \item{line.col}{the line color of the fitted distribution}
  \item{col.main}{the colour for the main title}
  \item{col.lab}{the colour of the labels}
  \item{col.axis}{the color of the axis}
 \item{\dots}{for extra arguments to be passed to the \code{gamlss} function}
}
\details{
  This function first fits constants for each parameters of a GAMLSS distribution family using the \code{gamlss} function  
  and them plots the fitted distribution together with the appropriate plot according to whether 
  the \code{y} variable is of a continuous or discrete type. Histogram is plotted for continuous and barplot for discrete variables. 
  The function \code{truehist} of 
  Venables and Ripley's MASS package is used for the histogram plotting.  
}
\value{
  returns a plot
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

\author{Mikis Stasinopoulos}

\seealso{\code{\link{gamlss}}, \code{\link[gamlss.dist]{gamlss.family}} }
\examples{
data(abdom)
histDist(y,family="NO", data=abdom)
# use the ylim 
histDist(y,family="NO", ylim=c(0,0.005), data=abdom)
# bad fit use PE
histDist(y,family="PE",ymax=0.005, data=abdom, line.col="blue")
# discere data counts 
#  Hand at al. p150  Leptinotarsa decemlineata
y <- c(0,1,2,3,4,6,7,8,10,11)
freq <- c(33,12,5,6,5,2,2,2,1,2) 
histDist(y, "NBI", freq=freq)
# the same as 
histDist(rep(y,freq), "NBI")
}
\keyword{regression}% 
