\name{cs}

\alias{cs}
\alias{scs}
\alias{cs.control}

\title{Specify a Smoothing Cubic Spline Fit in a GAMLSS Formula}

\description{The functions \code{cs()} and \code{scs()} are using the cubic smoothing splines function \code{smooth.spline()} to do smoothing.  They take a vector and return it with several attributes. 
The vector is used in the construction of the model matrix. The functions do not do the smoothing, but assigns the attributes to the vector to aid gamlss in the smoothing.
The function doing the smoothing is \code{\link{gamlss.cs}()}.
This function use the \R{} function \code{smooth.spline()} which is then used by the backfitting function  \code{additive.fit()} which is based on the original GAM implementation described in Chambers and Hastie (1992).
The function \code{gamlss.scs()}  differs from the function \code{cs()} in that allows cross validation of the smoothing parameters unlike the \code{cs()} which fixes the effective degrees of freedom, \code{df}.  Note that the recommended smoothing function is now the function \code{pb()} which allows the estimation of the smoothing parameters using a local maximum likelihood. The function \code{pb()} is based on the penalised beta splines (P-splines) of Eilers and Marx (1996).
  
The (experimental) function \code{vc} is now defunct. For fitting varying coefficient models, Hastie and Tibshirani (1993) use the function \code{pvc()}. 

}
\usage{
cs(x, df = 3, spar = NULL, c.spar = NULL, control = cs.control(...), ...)
scs(x, df = NULL, spar = NULL, control = cs.control(...), ...)
cs.control(cv = FALSE, all.knots = TRUE, nknots = NULL, keep.data = TRUE,
               df.offset = 0, penalty = 1.4, control.spar = list(), ...)
}

\arguments{
  \item{x}{the univariate predictor, (or expression, that evaluates to a numeric vector). 
           For the function \code{vc} the x argument is the vector which has its (linear) coefficient change with \code{r}}
  \item{df}{the desired equivalent number of degrees of freedom (trace of the smoother matrix minus two for the constant and linear fit). The real smoothing parameter (spar below) is found such that df=tr(S)-2, 
  where S is the implicit smoother matrix. Values for df should be greater than 0, with 0 implying a linear fit.}
  \item{spar}{ smoothing parameter, typically (but not necessarily) in (0,1].   
               The coefficient lambda of the integral of the squared second derivative in the fit (penalised log likelihood)
               criterion is a monotone function of `spar', see the details in \code{smooth.spline}.}
  \item{c.spar}{This is an option to be used when the degrees of freedom of the fitted gamlss object are different from the ones given as input in the option \code{df}.  The default values used are the ones given  the  option \code{control.spar}    in the R function \code{smooth.spine()} and they are \code{c.spar=c(-1.5, 2)}. For  very large data sets e.g. 10000 observations, the upper limit may have to increase for example to  \code{c.spar=c(-1.5, 2.5)}. Use this option if you have received the warning 'The output df are different from the input, change the control.spar'. \code{c.spar} can take both vectors or lists of length 2, for example \code{c.spar=c(-1.5, 2.5)} or \code{c.spar=list(-1.5, 2.5)} would have the same effect.}
  \item{control}{control for the function \code{smooth.spline()}, see below}
  \item{cv}{see the  \R{} function \code{smooth.spline()}}
  \item{all.knots}{see the  \R{} function \code{smooth.spline()}}
  \item{nknots}{see the  \R{} function \code{smooth.spline()}}
  \item{keep.data}{see the  \R{} function \code{smooth.spline()}}
  \item{df.offset}{see the  \R{} function \code{smooth.spline()}}
  \item{penalty}{see the  \R{} function \code{smooth.spline()}, here the default value is 1.4}
  \item{control.spar}{see above \code{c.spar} or the equivalent argument in the function \code{smooth.spline}}
  \item{\dots}{for extra arguments}
}
\details{
Note that \code{cs} itself does no smoothing; it simply sets things up for the function \code{gamlss()} which in turn uses the function 
  \code{additive.fit()} for backfitting which in turn uses \code{gamlss.cs()}

Note that  \code{cs()} and  \code{scs()} functions behave differently at their default values that is if df and lambda are not specified.
\code{cs(x)} by default will use 3 extra degrees of freedom for smoothing for \code{x}.
\code{scs(x)} by default will estimate lambda (and the degrees of freedom) automatically using generalised cross validation (GCV).
Note that if GCV is used the convergence of the gamlss model  can be less stable compared to a model where the degrees of freedom are fixed. This will be true for small data sets.         
  }


\value{
  the vector x is returned, endowed with a number of attributes. The vector itself is used in the construction of the model matrix, 
  while the attributes are needed for the backfitting algorithms \code{additive.fit()}. 
  Since smoothing splines includes linear fits, the linear part will be efficiently computed with the other parametric linear parts of the model. 
}
\references{

Chambers, J. M. and Hastie, T. J. (1992) \emph{Statistical Models in S}, Wadsworth & Brooks/Cole.

Eilers, P. H. C. and Marx, B. D. (1996). Flexible smoothing with B-splines and penalties (with comments and rejoinder). \emph{Statist. Sci}, \bold{11}, 89-121.

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

\author{Mikis Stasinopoulos and Bob Rigby (see also the documentation of the function\code{smooth.spline()} for the  original authors of the cubic spline function.)}

\note{The degrees of freedom df are defined differently from that of the gam() function in S-plus. Here df are the additional degrees of freedom 
excluding the constant and the linear part of x. For example \code{df=4} in \code{gamlss()} is equivalent to \code{df=5} in \code{gam()} in S-plus }

\section{Warning}{For a user who wishes to compare the \code{gamlss()} results with the equivalent \code{gam()} 
results in S-plus: make sure when using S-plus that the convergence criteria epsilon and bf.epsilon in \code{control.gam()} are decreased sufficiently
to ensure proper convergence in S-plus. 
Also note that the degrees of freedom are defined on top of the linear term in \code{gamlss}, but on top of the constant term in S-plus, 
(so use an extra degrees of freedom in S-plus in order to obtain comparable results to those in \code{galmss}).

Change the upper limit of spar if you received the warning 'The output df are different from the input, change the control.spar'.

For large data sets do not use expressions, e.g. \code{cs(x^0.5)} inside the \code{gamlss} function command but evaluate the expression, 
e.g. nx=\eqn{x^0.5}{x^0.5}, first and then use \code{cs(nx)}.  
} 

\seealso{ \code{\link{gamlss}}, \code{\link{gamlss.cs}}, \code{\link{pb}}, \code{\link{pvc}}   }

\examples{
# cubic splines example
data(aids)
# fitting a smoothing cubic spline with 7 degrees of freedom
# plus the a quarterly  effect  
aids1<-gamlss(y~cs(x,df=7)+qrt,data=aids,family=PO) # 
aids2<-gamlss(y~scs(x,df=5)+qrt,data=aids,family=PO) # 
aids3<-gamlss(y~scs(x)+qrt,data=aids,family=PO) # using GCV 
with(aids, plot(x,y))
lines(aids$x,fitted(aids1), col="red")
lines(aids$x,fitted(aids3), col="green")
rm(aids1, aids2, aids3)
}

\keyword{regression}% 
