\name{.binom}
\alias{.binom}
\alias{.counts}
\alias{.gamlss.bi.list}
\alias{.gamlss.multin.list}
\alias{.gamlss.sm.list}
\alias{.real0to1}
\alias{.realAll}
\alias{.realline}
\alias{.realplus}

\docType{data}
\title{Lists used by GAMLSS}
\description{
Those lists are used in GAMLSS fits.
}
\usage{

".binom"

".counts"

".gamlss.bi.list" 

".gamlss.multin.list"

".gamlss.sm.lis"

".real0to1"

".realAll" 

".realline"

".realplus" 
}
\format{
  List used by the \code{gamlss()} function.
  \describe{
    \item{\code{.binom}}{a character vector showing all the binomial type (finite count) distributions}
     \item{\code{.counts}}{a character vector showing all the infinity count distributions}
    \item{\code{.gamlss.bi.list}}{a character vector showing all the binomial type (finite count) distributions}
    \item{\code{.gamlss.multin.list}}{a character vector showing all the multinomial distributions}
    \item{\code{.gamlss.sm.list}}{a character vector showing all the available smooth functions}
    \item{\code{.real0to1}}{a character vector showing all real line distributions with range 0 to 1 }
    \item{\code{.realAll}}{a character vector showing all real line distributions from \eqn{-\infty} to \eqn{+\infty} and from \eqn{0} to \eqn{+\infty}}
    \item{\code{.realline}}{a character vector showing all all real line distributions from \eqn{-\infty} to \eqn{+\infty}}
    \item{\code{.realplus}}{a character vector showing all all real line distributions from \eqn{0} to \eqn{+\infty}}
  }
}
\details{Those list are internal to help the \code{gamlss()} function. 
}

\references{
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
\examples{
.binom
.counts
.gamlss.bi.list
.gamlss.multin.list
.gamlss.sm.list
.real0to1
.realline
.realplus
}
\keyword{datasets}
