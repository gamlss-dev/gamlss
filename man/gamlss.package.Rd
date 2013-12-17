\name{gamlss-package}
\alias{gamlss-package}
\docType{package}
\title{
The GAMLSS library and datasets
}
\description{
This a collection of functions to fit Generalized Additive Models for Location Scale and Shape(GAMLSS)
and handled gamlss objects. 

GAMLSS were introduced by Rigby and Stasinopoulos (2005).
GAMLSS is a general framework for univariate regression type statistical
problems using new ways of dealing with overdispersion, skewness
and kurtosis in the response variable.
In GAMLSS the exponential family distribution assumption used in
Generalized Linear Model (GLM) and Generalized Additive Model
(GAM),(see Nelder and Wedderburn, 1972 and Hastie and Tibshirani,
1990, respectively) is relaxed and replaced by a very general
distribution family including highly skew and kurtotic discrete
and continuous distributions. The systematic part of the model is
expanded to allow modelling not only the mean (or location) but
other parameters of the distribution of the response variable as linear parametric,
nonlinear parametric or additive non-parametric functions of
explanatory variables and/or random effects terms. Maximum
(penalized) likelihood estimation is used to fit the models. 
}
\details{
\tabular{ll}{
Package: \tab gamlss\cr
Type: \tab Package\cr
Version: \tab 1.5-0\cr
Date: \tab 2006-12-13\cr
License: \tab GPL (version 2 or later) See file LICENSE\cr
}
This package allow the user to model the distribution of the response variable 
using a variety of one, two, three and four parameter families of
distributions. The distributions implemented currently can be found in \code{\link[gamlss.dist]{gamlss.family}}.
Other distributions can be easily added.  
In the current implementation of GAMLSS several additive terms
have been implemented including regression splines, smoothing
splines, penalized splines, varying coefficients, fractional
polynomials and random effects. Other additive terms can be easily
added.  
}
\author{Mikis Stasinopoulos <d.stasinopoulos@londonmet.ac.uk>, Bob Rigby <r.rigby@londonmet.ac.uk> with
        contributions from Calliope Akantziliotou.

Maintainer: Mikis Stasinopoulos <d.stasinopoulos@londonmet.ac.uk>
}
\references{

Nelder, J. A. and Wedderburn, R. W. M. (1972). Generalized
linear models. \emph{J. R. Statist. Soc.} A., \bold{135} 370-384.

Hastie, T. J. and Tibshirani, R. J. (1990). \emph{Generalized Additive
Models}. Chapman and Hall, London.

Rigby, R. A. and  Stasinopoulos D. M. (2005). Generalized additive models for location, scale and shape,(with discussion), 
\emph{Appl. Statist.}, \bold{54}, part 3, pp 507-554.

Stasinopoulos D. M., Rigby R.A. and Akantziliotou C. (2006) Instructions on how to use the GAMLSS package in R.
Accompanying documentation in the current GAMLSS  help files, (see also  \url{http://www.gamlss.org/}).

Stasinopoulos D. M. Rigby R.A. (2007) Generalized additive models for location scale and shape (GAMLSS) in R.
\emph{Journal of Statistical Software}, Vol. \bold{23}, Issue 7, Dec 2007, \url{http://www.jstatsoft.org/v23/i07}.
}
\keyword{ package }

\examples{
data(abdom)
mod<-gamlss(y~pb(x),sigma.fo=~pb(x),family=BCT, data=abdom, method=mixed(1,20))
plot(mod)
rm(mod)
}
