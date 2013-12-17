\name{bfp}
\alias{bfp}
\alias{fp}
\alias{pp}


\title{ Functions to fit fractional polynomials in GAMLSS}
\description{
The function \code{bfp} generate a power polynomial basis matrix which  (for given powers) can be used to fit power polynomials in one x-variable.
The function \code{fp} takes a vector and returns it with several attributes. 
The vector is used in the construction of the model matrix.  The function \code{fp()} is not used for fitting the fractional polynomial curves
but assigns the attributes to the vector to aid gamlss in the fitting process.
The function doing the fitting is \code{\link{gamlss.fp}()} which is used at the backfitting function 
\code{\link{additive.fit}} (but never used on its own).
The (experimental) function \code{pp} can be use to fit power polynomials as in \eqn{a+b_1 
x^{p_1}+b_2 x^{p_2}}{a+b1*x^p1+b2*x^p2}., where p1 and p2 
have arbitrary values rather restricted as in the  \code{fp} function.   

}
\usage{
bfp(x, powers = c(1, 2), shift = NULL, scale = NULL)
fp(x, npoly = 2, shift = NULL, scale = NULL)
pp(x, start = list(), shift = NULL, scale = NULL)
}

\arguments{
  \item{x}{the explanatory variable to be used in functions \code{bfp()} or \code{fp()}. Note that this is different from the argument \code{x} use in \code{gamlss.fp} (a function used in the backfitting but not by straight by the user) }
  \item{powers}{a vector containing as elements  the powers in which the x has to be raised }
  \item{shift}{a number for shifting the x-variable. The default values is  zero, if x is positive, or the minimum of the positive difference in x  minus the minimum of x}
  \item{scale}{a positive number for scalling the x-variable. The default values is 
               \eqn{10^(sign(log10(range)))*trunc(abs(log10(range)))}{10^(sign(log10(range)))*trunc(abs(log10(range)))}}
  \item{npoly}{a positive indicating how many fractional polynomials should be considered in the fit. Can take the values 1, 2 or 3 with 2 as default}
  \item{start}{a list containing the starting values for the non-linear maximization to find the powers. 
             The results from fitting the equivalent fractional polynomials can be used here  }
 } 

\details{ The above functions are an implementation of the 
fractional polynomials introduced by  Royston and Altman (1994). 
The three functions involved in the fitting are loosely based on 
the fractional polynomials implementation in S-plus written by 
Gareth Amber.  The function \code{bfp} generates the right design 
matrix for the fitting a power polynomial of the type \eqn{a+b_1 
x^{p_1}+b_2 x^{p_2}+\ldots+b_k x^p_k }{a+b1*x^p1+b2*x^p2+...+bk*x^pk}. For given powers 
\eqn{p_1,p_2,\ldots,p_k}{p1,p2,...,pk} given as the argument \code{powers} 
in \code{bfp()} the function can be used to fit power polynomials 
in the same way as the functions \code{poly()} or \code{bs()} (of 
package \code{splines}) are used to fit orthogonal or piecewise 
polynomials respectively.    
The function \code{fp()}, which is working as a smoother in \code{gamlss}, is used to fit the best fractional polynomials within a set of power values. 
Its argument \code{npoly} 
determines whether one, two or three fractional polynomials should used in the fitting. 
For a fixed number \code{npoly} the algorithm looks for the best fitting fractional polynomials
in the list \code{c(-2, -1, -0.5, 0, 0.5, 1, 2, 3) }. Note that \code{npolu=3} is rather slow since it fits all possible combinations 3-way combinations 
at each backfitting interaction. 
The function \code{gamlss.fp()} is an internal function of GAMLSS allowing the 
fractional polynomials to be fitted in the backfitting cycle of \code{gamlss}, and should be not used on its own.  
}
\value{
  The function \code{bfp} returns a matrix to be used as part of the design matrix in the fitting. 
  
  The function \code{fp} returns a vector with values zero to be included in the design matrix but with attributes useful in the fitting 
  of the fractional polynomials algorithm in \code{gamlss.fp}.  

}
\references{

Amber G. (1999) Fracial polynomials in S-plus,  \url{http://lib.stat.cmu.edu/S/fracpoly}.

Rigby, R. A. and  Stasinopoulos D. M. (2005). Generalized additive models for location, scale and shape,(with discussion), 
\emph{Appl. Statist.}, \bold{54}, part 3, pp 507-554.

Royston,  P. and Altman, D. G., (1994). Regression using fractional polynomials of continuous covariates: parsimonious parametric modelling 
                   (with discussion), \emph{Appl. Statist.}, \bold{43}, 429-467.
                   
Stasinopoulos D. M., Rigby R.A. and Akantziliotou C. (2006) Instructions on how to use the GAMLSS package in R.
Accompanying documentation in the current GAMLSS  help files, (see also  \url{http://www.gamlss.org/}).

Stasinopoulos D. M. Rigby R.A. (2007) Generalized additive models for location scale and shape (GAMLSS) in R.
\emph{Journal of Statistical Software}, Vol. \bold{23}, Issue 7, Dec 2007, \url{http://www.jstatsoft.org/v23/i07}.
}
 
\author{Mikis Stasinopoulos \email{d.stasinopoulos@londonmet.ac.uk}, Bob Rigby \email{b.rigby@londonmet.ac.uk} }


\section{Warning }{Since the model constant is included in both the design matrix X and in the backfitting part of fractional polynomials, its values is wrongly 
  given in the \code{summary}. Its true values is the model constant minus the constant from the fractional polynomial fitting ??? What happens if more that one fractional polynomials are fitted? } 

\seealso{ \code{\link{gamlss}}, \code{\link[gamlss.dist]{gamlss.family}} }
\examples{
data(abdom)
#fits polynomials with power 1 and .5 
mod1<-gamlss(y~bfp(x,c(1,0.5)),data=abdom)
# fit the best of one fractional polynomial
m1<-gamlss(y~fp(x,1),data=abdom)
# fit the best of two fractional polynomials
m2<-gamlss(y~fp(x,2),data=abdom)
# fit the best of three fractional polynomials
m3<-gamlss(y~fp(x,3),data=abdom)
# get the coefficient for the second model 
m2$mu.coefSmo
# now power polynomials using the best 2 fp c()
 m4 <- gamlss(y ~ pp(x, c(1,3)), data = abdom)
# This is not good idea in this case because
# if you look at the fitted values you see what it went wrong
plot(y~x,data=abdom)
lines(fitted(m2,"mu")~abdom$x,col="red")
lines(fitted(m4,"mu")~abdom$x,col="blue")
}

\keyword{regression}%
