\name{centiles.split}
\alias{centiles.split}
\title{Plots centile curves split by x for a GAMLSS object}
\description{
This function plots centiles curves for separate ranges of the unique explanatory variable x. 
It is similar to the \code{\link{centiles}} function but the range of x is split at a user defined values \code{xcut.point} into r separate ranges.
The functions also tabulates the sample percentages below each centile curve for each of the r ranges of x 
(for comparison with the model percentage given by cent)
 The model should have only one explanatory variable.    
}

\usage{
centiles.split(obj, xvar = NULL, xcut.points = NULL, n.inter = 4, 
               cent = c(0.4, 2, 10, 25, 50, 75, 90, 98, 99.6), 
               legend = FALSE, main = NULL, main.gsub = "@", 
               ylab = "y", xlab = "x", ylim = NULL, overlap = 0, 
               save = TRUE, plot = TRUE, ...)
}

\arguments{
  \item{obj}{a fitted gamlss object from fitting a gamlss continuous distribution}
  \item{xvar}{the unique explanatory variable}
  \item{xcut.points}{the x-axis cut off points e.g. \code{c(20,30)}. If \code{xcut.points=NULL} then the \code{n.inter} argument is activated }
  \item{n.inter}{if \code{xcut.points=NULL} this argument gives the number of intervals in which the x-variable will be splited, with default 4 }
  \item{cent}{a vector with elements the \% centile values for which the centile curves are to be evaluated}
  \item{legend}{whether a legend is required in the plots or not, the default is \code{legent=FALSE}  }
  \item{main}{the main title as character. 
              If NULL the default title (shown the intervals) is shown}
  \item{main.gsub}{if the \code{main.gsub} (with default "@") appears in the \code{main} 
                 title then it is substituted with the default title.}
  \item{ylab}{the y-variable label}
  \item{xlab}{the x-variable label}
  \item{ylim}{the range of the y-variable axis}
  \item{overlap}{how much overlapping in the \code{xvar} intervals. Default value is \code{overlap=0} for non overlapping intervals}
  \item{save}{whether to save the sample percentages or not with default equal to \code{TRUE}. 
            In this case the functions produce a matrix giving the sample percentages for each interval }
  \item{plot}{whether to plot the centles. This option is usefull if the sample statistics only are to be used}
  \item{\dots}{for extra arguments}
}

 
\value{Centile plots are produced and the sample centiles below each centile curve for each of the r ranges of x can be saved into a matrix.  
}
\references{
Rigby, R. A. and  Stasinopoulos D. M. (2005). Generalized additive models for location, scale and shape,(with discussion), 
\emph{Appl. Statist.}, \bold{54}, part 3, pp 507-554.

Stasinopoulos D. M., Rigby R.A. and Akantziliotou C. (2003) Instructions on how to use the GAMLSS package in R.
Accompanying documentation in the current GAMLSS  help files, (see also  \url{http://www.gamlss.org/}).

Stasinopoulos D. M. Rigby R.A. (2007) Generalized additive models for location scale and shape (GAMLSS) in R.
\emph{Journal of Statistical Software}, Vol. \bold{23}, Issue 7, Dec 2007, \url{http://www.jstatsoft.org/v23/i07}.
}
\author{Mikis Stasinopoulos, \email{d.stasinopoulos@londonmet.ac.uk}, Bob Rigby \email{r.rigby@londonmet.ac.uk}, with contributions from Elaine Borghie}


\section{Warning }{This function is appropriate when only one continuous explanatory variable is fitted in the model}
\seealso{\code{\link{gamlss}}  \code{\link{centiles}},   \code{\link{centiles.com}} }
\examples{
data(abdom)
h<-gamlss(y~pb(x), sigma.formula=~pb(x), family=BCT, data=abdom) 
mout <- centiles.split(h,xvar=abdom$x)
mout
rm(h,mout)
}
\keyword{regression}% 
