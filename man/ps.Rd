\name{ps}
\alias{ps}
\alias{pb}
\alias{pb.control}
\alias{pvc}
\alias{pvc.control}
\alias{cy}
\alias{cy.control}

\title{Specify a Penalised  B-Spline Fit in a GAMLSS Formula}
\description{There are several function operating using penalised B-splines: \code{pb()},  \code{cy()}, \code{pvc()} and \code{ps()}.
The functions take a vector and return it with several attributes. 
The vector is used in the construction of the design matrix X used in the fitting.  The functions  do not do the  smoothing, 
but assign the attributes to the vector to aid gamlss in the smoothing. 
The functions doing the smoothing are  \code{\link{gamlss.ps}()} \code{\link{gamlss.pb}()}, \code{\link{gamlss.cy}()} and \code{\link{gamlss.pvc}()}
which are used in the backfitting  function \code{\link{additive.fit}}.

The function \code{pb()} is more efficient and faster than the original penalized smoothing function \code{ps()}. 
\code{pb()} allows the estimation of the smoothing parameters using different local (performance iterations) methods.
The method are "ML", "ML-1", "EM", "GAIC" and "GCV". The function \code{cy()} fits a cycle penalised beta regression spline such as 
the last fitted value of the smoother is equal to the first fitted value.
The function \code{pvc()} fits  varying coefficient models see Hastie and Tibshirani(1993) and it is more general and flexible than the old 
\code{vc()} function which is based on cubic splines.
}
\usage{
pb(x, df = NULL, lambda = NULL, control = pb.control(...), ...)
pb.control(inter = 20, degree = 3, order = 2, start = 10, quantiles = FALSE, 
               method = c("ML", "GAIC", "GCV", "EM", "ML-1"), k = 2, ...)
cy(x, df = NULL, lambda = NULL, control = cy.control(...), ...)
cy.control(inter = 20, degree = 3, order = 2, start = 10, 
          method = c("ML", "GAIC", "GCV", "EM", "ML-1"), k = 2, ts=FALSE, ...)
pvc(x, df = NULL, lambda = NULL, by = NULL, control = pvc.control(...), ...)          
pvc.control(inter = 20, degree = 3, order = 2, start = 10, quantiles = FALSE, 
             method = c("ML", "GAIC", "GCV", "EM", "ML-1"), k = 2, ...)          
ps(x, df = 3, lambda = NULL, ps.intervals = 20, degree = 3, order = 3)
}


\arguments{
  \item{x}{the univariate predictor}
  \item{df}{the desired equivalent number of degrees of freedom (trace of the smoother matrix minus two for the constant and linear fit)}
  \item{lambda}{the smoothing parameter}
  \item{control}{setting the control parameters}
  \item{by}{a factor, for fitting different smoothing curves to each level of the factor or a continuous explanatory variable in which case  
             the coefficients of the \code{by} variable change smoothly according to \code{x} i.e. beta(x)*z where z is the \code{by} variable. }
  \item{\dots}{for extra arguments}
  \item{inter}{the no of break points (knots) in the x-axis}
  \item{degree}{the degree of the piecewise polynomial}
  \item{order}{the required difference in the vector of coefficients}
  \item{start}{the lambda starting value if the local methods are used, see below}
  \item{quantiles}{if TRUE the quantile values of x are use to determine the knots} 
  \item{ts}{if TRUE assumes that it is a seasonal factor} 
  \item{method}{The method used in the (local) performance iterations. Available methods are "ML", "ML-1", "EM", "GAIC" and "GCV"}
  \item{k}{the penalty used in "GAIC" and "GCV"}
  \item{ps.intervals}{the no of break points in the x-axis}
} 
\details{
The \code{ps()} function is based on Brian Marx function which can be found in \url{http://www.stat.lsu.edu/faculty/marx/}.
The \code{pb()}, \code{cy()} and \code{pvc()} functions are based on Paul Eilers's original R functions. 
Note that  \code{ps()} and  \code{pb()} functions behave differently at their default values if df and lambda are not specified.
\code{ps(x)} by default  uses 3 extra degrees of freedom for smoothing \code{x}.
\code{pb(x)} by default  estimates lambda (and therefore the degrees of freedom) automatically using a "local" method.
The local (or performance iterations) methods available are: 
(i) local Maximum Likelihood, "ML", 
(ii) local Generalized Akaike information criterion, "GAIC",
(iii) local Generalized Cross validation "GCV" 
(iv) local EM-algorithm, "EM" (which is very slow) and 
(v) a modified version of the ML, "ML-1" which produce identical results with "EM" but faster.
  
Note that the local (or performance iterations) methods can occasionally  make the convergence of gamlss less stable
 compared to  models where the degrees of freedom are fixed.           
}

\value{
 the vector x is returned, endowed with a number of attributes. The vector itself is used in the construction of the model matrix, 
  while the attributes are needed for the backfitting algorithms \code{additive.fit()}.  
}
\references{

\url{http://www.stat.lsu.edu/faculty/marx/} 

Eilers, P. H. C. and Marx, B. D. (1996). Flexible smoothing with
B-splines and penalties (with comments and rejoinder). \emph{Statist. Sci},
\bold{11}, 89-121.

Hastie, T. J. and Tibshirani, R. J. (1993), Varying coefficient models (with discussion),J. R. Statist. Soc. B., \bold{55},
    757-796.


Rigby, R. A. and  Stasinopoulos D. M. (2005). Generalized additive models for location, scale and shape,(with discussion), 
\emph{Appl. Statist.}, \bold{54}, part 3, pp 507-554.

Stasinopoulos D. M. Rigby R.A. (2007) Generalized additive models for location scale and shape (GAMLSS) in R.
\emph{Journal of Statistical Software}, Vol. \bold{23}, Issue 7, Dec 2007, \url{http://www.jstatsoft.org/v23/i07}.
}
\author{ Mikis Stasinopoulos \email{d.stasinopoulos@londonmet.ac.uk}, Bob Rigby \email{r.rigby@londonmet.ac.uk } and Paul Eilers}

\section{Warning}{There are occations where the automatic local methods do not work. One accation which came to our attention is  when 
the range of the response variable values is very large. Scalling the response variable will solve the problem.} 

\seealso{ \code{\link{gamlss}}, \code{\link{gamlss.ps}}, \code{\link{cs}}}
\examples{
#==============================
# pb() and ps() functions
data(aids)
# fitting a smoothing cubic spline with 7 degrees of freedom
# plus the a quarterly  effect  
aids1<-gamlss(y~ps(x,df=7)+qrt,data=aids,family=PO) # fix df's 
aids2<-gamlss(y~pb(x,df=7)+qrt,data=aids,family=PO) # fix df's
aids3<-gamlss(y~pb(x)+qrt,data=aids,family=PO) # estimate lambda
with(aids, plot(x,y))
with(aids, lines(x,fitted(aids1),col="red"))
with(aids, lines(x,fitted(aids2),col="green"))
with(aids, lines(x,fitted(aids1),col="yellow"))
rm(aids1, aids2, aids3)
#=============================
\dontrun{
# cy()
# simulate data
set.seed(555)
x = seq(0, 1, length = 100)
y = sign(cos(1 * x * 2 * pi + pi / 4)) + rnorm(length(x)) * 0.2
plot(y~x)
m1<-gamlss(y~cy(x)) 
lines(fitted(m1)~x)
rm(y,x,m1)
#=============================
# the pvc() function
# function to generate data
genData <- function(n=200)
 {
f1 <- function(x)-60+15*x-0.10*x^2
f2 <- function(x)-120+10*x+0.08*x^2
set.seed(1441)
x1 <- runif(n/2, min=0, max=55)
x2 <- runif(n/2, min=0, max=55)
y1 <- f1(x1)+rNO(n=n/2,mu=0,sigma=20)
y2 <- f2(x2)+rNO(n=n/2,mu=0,sigma=30)
 y <- c(y1,y2)
 x <- c(x1,x2)
 f <- gl(2,n/2)
da<-data.frame(y,x,f)
da
}
da<-genData(500)
plot(y~x, data=da, pch=21,bg=c("gray","yellow3")[unclass(f)])
# fitting models
# smoothing x
m1 <- gamlss(y~pb(x), data=da)
# parallel smoothing lines
m2 <- gamlss(y~pb(x)+f, data=da)
# linear interaction
m3 <- gamlss(y~pb(x)+f*x, data=da)
# varying coefficient model
m4 <- gamlss(y~pvc(x, by=f), data=da)
GAIC(m1,m2,m3,m4)
# plotting the fit
lines(fitted(m4)[da$f==1][order(da$x[da$f==1])]~da$x[da$f==1]
         [order(da$x[da$f==1])], col="blue", lwd=2)
lines(fitted(m4)[da$f==2][order(da$x[da$f==2])]~da$x[da$f==2]
         [order(da$x[da$f==2])], col="red", lwd=2)
rm(da,m1,m2,m3,m4)
#=================================
# the rent data
# first with a factor
data(rent)
plot(R~Fl, data=rent, pch=21,bg=c("gray","blue")[unclass(rent$B)])
r1 <- gamlss(R~pb(Fl), data=rent)
# identical to model
r11 <- gamlss(R~pvc(Fl), data=rent)
# now with the factor
r2 <- gamlss(R~pvc(Fl, by=B), data=rent)
lines(fitted(r2)[rent$B==1][order(rent$Fl[rent$B==1])]~rent$Fl[rent$B==1]
                [order(rent$Fl[rent$B==1])], col="blue", lwd=2)
lines(fitted(r2)[rent$B==0][order(rent$Fl[rent$B==0])]~rent$Fl[rent$B==0]
                [order(rent$Fl[rent$B==0])], col="red", lwd=2)
# probably not very sensible model
rm(r1,r11,r2)
#-----------
# now with a continuous variable
# additive model
 h1 <-gamlss(R~pb(Fl)+pb(A), data=rent)
# varying-coefficient model
 h2 <-gamlss(R~pb(Fl)+pb(A)+pvc(A,by=Fl), data=rent)
AIC(h1,h2)
rm(h1,h2)}
#==================================
}
\keyword{regression}% 
