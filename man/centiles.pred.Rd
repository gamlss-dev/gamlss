\name{centiles.pred}
\alias{centiles.pred}

\title{Creating  predictive centiles values}
\description{
 This function creates predictive centiles curves for new x-values given a GAMLSS fitted model.
  The function has three options: i) for given new x-values and given percentage centiles calculates a matrix containing 
  the centiles values for y,
  ii) for given new x-values and standard normalized centile values calculates a matrix containing the centiles values for y,
  iii) for given new x-values and new y-values calculates the z-scores.  
A restriction of the function is that it applies to models with only one explanatory variable. 
}
\usage{
centiles.pred(obj, type = c("centiles", "z-scores", "standard-centiles"), 
             xname = NULL, xvalues = NULL, power = NULL, yval = NULL, 
             cent = c(0.4, 2, 10, 25, 50, 75, 90, 98, 99.6), 
             dev = c(-4, -3, -2, -1, 0, 1, 2, 3, 4), calibration = FALSE,
             plot = FALSE, legend = TRUE,  ylim = NULL,xlim = NULL,
             ...)
}

\arguments{
  \item{obj}{a fitted gamlss object from fitting a gamlss continuous distribution }
  \item{type}{the default, "centiles", gets the centiles values given in the option \code{cent}. 
             \code{type="standard-centiles"} gets the standard centiles  given in the \code{dev}. 
             \code{type="z-scores"} gets the z-scores for given y and x new values}
  \item{xname}{ the name of the unique explanatory variable (it has to be the same as in the original fitted model)}
  \item{xvalues}{the new values for the explanatory variable where the prediction will take place}
  \item{power}{if power transformation is needed (but read the note below)}
  \item{yval}{the response values for a given x required for the calculation of "z-scores"}
  \item{cent}{a vector with elements the \% centile values for which the centile curves have to be evaluated}
  \item{dev}{ a vector with elements the standard normalized values for which the centile curves have to be evaluated in the option \code{type="standard-centiles"}}
   \item{calibration}{whether to calibrate the "centiles", the default is \code{calibrate=FALSE}}
  \item{plot}{whether to plot the "centiles" or the "standard-centiles", the default is \code{plot=FALSE}}
  \item{legend}{whether a legend is required in the plot or not, the default is \code{legent=TRUE}  } 
  \item{ylim}{If different \code{ylim} is required from the default} 
  \item{xlim}{If different \code{xlim} is required from the default } 
  \item{\dots}{ for extra arguments }
}

\value{
  a vector (for option \code{type="z-scores"}) or a  matrix for options
   \code{type="centiles"} or \code{type="standard-centiles"}
  containing the appropriate values
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
\author{Mikis Stasinopoulos, based on ideas of Elaine Borghie 
from the World Health Organization}
\note{The power option should be only used if the model  }

\section{Warning }{See example below of how to use the function when power transformation is used for the x-variables} 

\seealso{\code{\link{gamlss}}, \code{\link{centiles}}, \code{\link{centiles.split}} }
\examples{
## bring the data and fit the model
data(abdom)
a<-gamlss(y~pb(x),sigma.fo=~pb(x), data=abdom, family=BCT)
## plot the centiles
centiles(a,xvar=abdom$x)
##-----------------------------------------------------------------------------
## the first use of the function centiles.pred()
## to calculate the centiles at new x values
##-----------------------------------------------------------------------------
newx<-seq(12,40,2)
mat <- centiles.pred(a, xname="x", xvalues=newx )
mat
## now plot the centile curves  
 mat <- centiles.pred(a, xname="x",xvalues=newx, plot=TRUE )
##-----------------------------------------------------------------------------
## the second use of the function centiles.pred()
## to calculate (nornalised) standard-centiles for new x
## values using the fitted model
##-----------------------------------------------------------------------------
newx <- seq(12,40,2)
mat <- centiles.pred(a, xname="x",xvalues=newx, type="standard-centiles" )
mat
## now plot the standard centiles  
mat <- centiles.pred(a, xname="x",xvalues=newx, type="standard-centiles",
       plot = TRUE )
##-----------------------------------------------------------------------------
## the third use of the function centiles.pred()
##  if we have new x and y values what are their z-scores?
##-----------------------------------------------------------------------------
# create new y and x values and plot them in the previous plot
newx <- c(20,21.2,23,20.9,24.2,24.1,25)
newy <- c(130,121,123,125,140,145,150)
for(i in 1:7) points(newx[i],newy[i],col="blue")
## now calculate their z-scores
znewx <- centiles.pred(a, xname="x",xvalues=newx,yval=newy, type="z-scores" )
znewx
\dontrun{
##-----------------------------------------------------------------------------
## What we do if the x variables is transformed?
##----------------------------------------------------------------------------
##  case 1 : transformed x-variable within the formula
##----------------------------------------------------------------------------
## fit model
aa <- gamlss(y~pb(x^0.5),sigma.fo=~pb(x^0.5), data=abdom, family=BCT)
## centiles is working in this case
centiles(aa, xvar=abdom$x, legend = FALSE)
## get predict for values of x at 12, 14, ..., 40
mat <- centiles.pred(aa, xname="x", xvalues=seq(12,40,2), plot=TRUE )
mat
# plot all prediction points
xx <- rep(mat[,1],9)
yy <- unlist(mat[,2:10])
points(xx,yy,col="red")
##----------------------------------------------------------------------------
##  case 2 : the x-variable is previously transformed 
##----------------------------------------------------------------------------
nx <- abdom$x^0.5
aa <- gamlss(y~pb(nx),sigma.fo=~pb(nx), data=abdom, family=BCT)
centiles(aa, xvar=abdom$x)
# equivalent to fitting
newd<-data.frame( abdom, nx=abdom$x^0.5)
aa1 <- gamlss(y~pb(nx),sigma.fo=~pb(nx), family=BCT, data=newd)
centiles(aa1, xvar=abdom$x)
# getting the centiles at x equal to 12, 14, ...40
mat <-  centiles.pred(aa, xname="nx", xvalues=seq(12,40,2), power=0.5, 
         data=newd, plot=TRUE)
# plot all prediction points         
xxx <- rep(mat[,1],9)
yyy <- unlist(mat[,2:10])
points(xxx,yyy,col="red")
# the idea is that if the transformed x-variable is used in the fit
# the power argument has to used in centiles.pred()
}
}
\keyword{regression}% 
