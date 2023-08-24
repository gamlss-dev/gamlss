\name{ps}
\alias{ps}
\alias{pb}
\alias{pbo}
\alias{pbp}
\alias{pbc}
\alias{pbc.control}
\alias{pb.control}
\alias{pbo.control}
\alias{pbp.control}
\alias{pvc}
\alias{pvc.control}
\alias{cy}
\alias{cy.control}
\alias{pbm}
\alias{pbm.control}
\alias{pbz}
\alias{pbz.control}
\alias{getZmatrix}


\title{P-Splines Fits in a GAMLSS Formula}

\description{
There are several function which use P-spline methodology:

a) \code{pb()}, the current version of P-splines which uses SVD in the fitting and therefore is the most reliable 

b) \code{pbo()} and \code{pbp()}, older versions of P-splines. The first  uses a simple matrix algebra in the fits. The second is the last version of \code{pb()} with SVD  but uses different method  for prediction.

c) \code{pbc()} the new version of cycle P-splines (using SVD)

d) \code{cy()} the older version of cycle P-splines.

e) \code{pbm()} for fitting monotonic  P-splines (using SVD)

f) \code{pbz()} for fitting  P-splines which allow the fitted curve to shrink to zero degrees of freedom 

g) \code{ps()} the original P-splines with no facility of estimating the smoothing parameters and 

j) \code{pvc()} penalised varying coefficient models.

k) \code{pvp()} older version of pb() where the prediction was different (it is here in case someone would like to compare the results).

Theoretical explanation of the above P-splines can be found in Eilers \emph{et al.} (2016)

The functions take a vector and return it with several attributes. The vector is used in the construction of the design matrix X used in the fitting. The functions do not do the  smoothing, but assign the attributes to the vector to aid gamlss in the smoothing. The functions doing the smoothing are   \code{\link{gamlss.pb}()}, \code{\link{gamlss.pbo}()}, \code{\link{gamlss.pbc}()}  \code{\link{gamlss.cy}()}  \code{\link{gamlss.pvc}()},  \code{\link{gamlss.pbm}()}, \code{\link{gamlss.pbz}} and \code{\link{gamlss.ps}()} which are used in the backfitting  function \code{\link{additive.fit}}.


The function \code{pb()} is more efficient and faster than the original penalised smoothing function \code{ps()}. After December 2014 the  \code{pb()} has changed radically to improved performance. The older version of the \code{pb()} function is called now \code{pbo()}. 
\code{pb()} allows the estimation of the smoothing parameters using different local (performance iterations) methods. The method are "ML", "ML-1", "EM", "GAIC" and "GCV". 

The function \code{pbm()} fits monotonic smooth functions, that is functions which increase or decrease monotonically depending on the value of the argument \code{mono} which takes the values \code{"up"} or \code{"down"}.

The function \code{pbz()} is similar to \code{pb()} with the extra property that when lambda becomes very large the resulting smooth function  goes to a constant rather than to a linear function. This is very useful for model selection. The function is based on Maria Durban idea of using a double penalty, one of order 2 and one of order 1. The second penalty only applies if the effective df are close to 2 (that is if a linear is already selected). 
  
The function \code{pbc()} fits a cycle penalised beta regression spline such as  the last fitted value of the smoother is equal to the first fitted value. \code{cy()} is the older version.

The function \code{pvc()} fits  varying coefficient models see Hastie and Tibshirani(1993) and it is more general and flexible than the old \code{vc()} function which was based on cubic splines.

The function \code{getZmatrix()} creates a (random effect) design matrix \code{Z} which can be  used to fit a P-splines smoother using the \code{re())} function. (The \code{re()}  is an interface with the random effect function \code{lme} of the package \pkg{nlme}. 
}
\usage{
pb(x, df = NULL, lambda = NULL, max.df=NULL, 
    control = pb.control(...), ...)
pbo(x, df = NULL, lambda = NULL, control = pbo.control(...), ...)
pbp(x, df = NULL, lambda = NULL, control = pbp.control(...), ...)
pbo.control(inter = 20, degree = 3, order = 2, start = 10, quantiles = FALSE, 
               method = c("ML", "GAIC", "GCV", "EM", "ML-1"), k = 2, ...)
pb.control(inter = 20, degree = 3, order = 2, start = 10, quantiles = FALSE, 
          method = c("ML", "GAIC", "GCV"), k = 2, ...)
pbp.control(inter = 20, degree = 3, order = 2, start = 10, quantiles = FALSE, 
          method = c("ML", "GAIC", "GCV"), k = 2, ...)
pbc(x,  df = NULL, lambda = NULL, max.df=NULL, 
    control = pbc.control(...), ...)
pbc.control(inter = 20, degree = 3, order = 2, start = 10, 
          method = c("ML", "GAIC", "GCV"), k = 2, sin = TRUE, ...)
cy(x, df = NULL, lambda = NULL, control = cy.control(...), ...)
cy.control(inter = 20, degree = 3, order = 2, start = 10, 
          method = c("ML", "GAIC", "GCV", "EM", "ML-1"), k = 2, ts=FALSE, ...)
pvc(x, df = NULL, lambda = NULL, by = NULL, control = pvc.control(...), ...)          
pvc.control(inter = 20, degree = 3, order = 2, start = 10, quantiles = FALSE, 
             method = c("ML", "GAIC", "GCV"), k = 2, ...) 
pbm(x, df = NULL, lambda = NULL, mono=c("up", "down"), 
            control = pbm.control(...), ...)
pbm.control(inter = 20, degree = 3, order = 2, start = 10, quantiles = FALSE, 
            method=c("ML","GAIC", "GCV"), k=2, kappa = 1e10, ...)
pbz(x, df = NULL, lambda = NULL, control = pbz.control(...), ...)
pbz.control(inter = 20, degree = 3, order = 2, start = c(1e-04, 1e-04), 
     quantiles = FALSE, method = c("ML", "GAIC", "GCV"), k = 2, lim = 3, ...)

ps(x, df = 3, lambda = NULL, ps.intervals = 20, degree = 3, order = 3)

getZmatrix(x, xmin = NULL, xmax = NULL, inter = 20, degree = 3, order = 2)
}

\arguments{
  \item{x}{the univariate predictor}
  \item{df}{the desired equivalent number of degrees of freedom (trace of the smoother matrix minus two for the constant and linear fit)}
  \item{lambda}{the smoothing parameter}
  \item{max.df}{the limit of how large the effective degrees of freedom should be allowed to be}
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
  \item{mono}{for monotonic P-splines whether going "up"  or "down"}
  \item{kappa}{the smoothing hyper-parameter for the monotonic part of smoothing}
  \item{ps.intervals}{the no of break points in the x-axis}
  \item{xmin}{minimum value for creating the B-spline}
  \item{xmax}{maximum value for creating the B-spline}
  \item{sin}{whether to use the sin penalty or not}
  \item{lim}{at which level the second penalty of order 1 should start}
} 
\details{
The \code{ps()} function is based on Brian Marx function which can be found in his website.
The \code{pb()}, \code{cy()}, \code{pvc()} and  \code{pbm()} functions are based on Paul Eilers's original R functions. 
Note that  \code{ps()} and  \code{pb()} functions behave differently at their default values if df and lambda are not specified.
\code{ps(x)} by default  uses 3 extra degrees of freedom for smoothing \code{x}.
\code{pb(x)} by default  estimates lambda (and therefore the degrees of freedom) automatically using a "local" method.
The local (or performance iterations) methods available are: 
(i) local Maximum Likelihood, "ML", 
(ii) local Generalized Akaike information criterion, "GAIC",
(iii) local Generalized Cross validation "GCV" 
(iv) local EM-algorithm, "EM" (which is very slow) and 
(v) a modified version of the ML, "ML-1" which produce identical results with "EM" but faster.

The  function \code{pb()} fits a P-spline smoother.

The  function \code{pbm()} fits a monotonic (going up or down) P-spline smoother.

The  function \code{pbc()} fits a P-spline smoother where the beginning and end are the same.  

The \code{pvc()} fits a varying coefficient model.  

Note that the local (or performance iterations) methods can occasionally  make the convergence of gamlss less stable compared to models where the degrees of freedom are fixed.           
}

\value{
 the vector x is returned, endowed with a number of attributes. The vector itself is used in the construction of the model matrix, 
  while the attributes are needed for the backfitting algorithms \code{additive.fit()}.  
}
\references{

Eilers, P. H. C. and Marx, B. D. (1996). Flexible smoothing with
B-splines and penalties (with comments and rejoinder). \emph{Statist. Sci},
\bold{11}, 89-121.

Eilers, Paul HC, Marx, Brian D and Durban, Maria, (2016) Twenty years of P-splines. \emph{SORT-Statistics and Operations Research Transactions}, \bold{39}, 149--186.


Hastie, T. J. and Tibshirani, R. J. (1993), Varying coefficient models (with discussion),J. R. Statist. Soc. B., \bold{55},
    757-796.


Rigby, R. A. and  Stasinopoulos D. M. (2005). Generalized additive models for location, scale and shape,(with discussion), 
\emph{Appl. Statist.}, \bold{54}, part 3, pp 507-554.

Rigby, R. A., Stasinopoulos, D. M.,  Heller, G. Z.,  and De Bastiani, F. (2019)
	\emph{Distributions for modeling location, scale, and shape: Using GAMLSS in R}, Chapman and Hall/CRC. An older version can be found in \url{https://www.gamlss.com/}.

Stasinopoulos D. M. Rigby R.A. (2007) Generalized additive models for location scale and shape (GAMLSS) in R.
\emph{Journal of Statistical Software}, Vol. \bold{23}, Issue 7, Dec 2007, \url{https://www.jstatsoft.org/v23/i07/}.

Stasinopoulos D. M., Rigby R.A., Heller G., Voudouris V., and De Bastiani F., (2017)
\emph{Flexible Regression and Smoothing: Using GAMLSS in R},  Chapman and Hall/CRC.  

(see also \url{https://www.gamlss.com/}).
}
\author{ Mikis Stasinopoulos, Bob Rigby and Paul Eilers}

\section{Warning}{There are occasions where the automatic local methods do not work. One accusation which came to our attention is  when 
the range of the response variable values is very large. Scaling the response variable will solve the problem.} 

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
# pbc()
# simulate data
set.seed(555)
x = seq(0, 1, length = 100)
y = sign(cos(1 * x * 2 * pi + pi / 4)) + rnorm(length(x)) * 0.2
plot(y~x)
m1<-gamlss(y~pbc(x)) 
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
rm(h1,h2)
#-----------
# monotone function
set.seed(1334)
x = seq(0, 1, length = 100)
p = 0.4
y = sin(2 * pi * p * x) + rnorm(100) * 0.1
plot(y~x)
m1 <- gamlss(y~pbm(x))
points(fitted(m1)~x, col="red")
yy <- -y
plot(yy~x)
m2 <- gamlss(yy~pbm(x, mono="down"))
points(fitted(m2)~x, col="red")
#==========================================
# the pbz() function
# creating uncorrelated data
set.seed(123)
y<-rNO(100)
x<-1:100
plot(y~x)
#----------------------
# ML estimation
m1<-gamlss(y~pbz(x))
m2 <-gamlss(y~pb(x))
AIC(m1,m2)
op <- par( mfrow=c(1,2))
term.plot(m1, partial=T)
term.plot(m2, partial=T)
par(op)
# GAIC estimation
m11<-gamlss(y~pbz(x, method="GAIC", k=2))
m21 <-gamlss(y~pb(x, method="GAIC", k=2))
AIC(m11,m21)
op <- par( mfrow=c(1,2))
term.plot(m11, partial=T)
term.plot(m21, partial=T)
par(op)
# GCV estimation
m12<-gamlss(y~pbz(x, method="GCV"))
m22 <-gamlss(y~pb(x, method="GCV"))
AIC(m12,m22)
op <- par( mfrow=c(1,2))
term.plot(m12, partial=T)
term.plot(m22, partial=T)
par(op)
# fixing df is more trycky since df are the extra df 
m13<-gamlss(y~pbz(x, df=0))
m23 <-gamlss(y~pb(x, df=0))
AIC(m13,m23)
# here the second penalty is not take effect therefore identical results 
m14<-gamlss(y~pbz(x, df=1))
m24 <-gamlss(y~pb(x, df=1))
AIC(m14,m24)
# fixing lambda
m15<-gamlss(y~pbz(x, lambda=1000))
m25 <-gamlss(y~pb(x, lambda=1000))
AIC(m15,m25)
#--------------------------------------------------
# prediction 
m1<-gamlss(y~pbz(x), data=data.frame(y,x))
m2 <-gamlss(y~pb(x), data=data.frame(y,x))
AIC(m1,m2)
predict(m1, newdata=data.frame(x=c(80, 90, 100, 110)))
predict(m2, newdata=data.frame(x=c(80, 90, 100, 110)))
#---------------------------------------------------
}
}
\keyword{regression}% 
