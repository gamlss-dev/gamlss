\name{gamlss.lo}
\alias{gamlss.lo}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Support for Function lo()}
\description{
 This is support for the \code{loess} function \code{lo()}.
It is not intended to be called directly by users. The function \code{gamlss.lo} is calling  the \R{} function \code{loess}.
}
\usage{
gamlss.lo(x, y, w, xeval = NULL, ...)
}
\arguments{
  \item{x}{the design matrix}
  \item{y}{the response variable}
  \item{w}{ prior weights}
   \item{xeval}{used in prediction } 
  \item{...}{further arguments passed to or from other methods.}
}

\value{

 Returns an object
 
  \item{fitted}{the smooth values}
  \item{residuals}{the residuals}
  \item{var}{the variance of the smoother}
  \item{nl.df}{the non-linear degrees of freedom}
  \item{coefSmo}{with value NULL}
  \item{lambda}{the value of span}
 
}

\author{Mikis Stasinopoulos based on Brian Ripley implementation of \code{loess} function in \R{}}


\seealso{ \code{\link{gamlss}}, \code{\link{lo}} }


\keyword{regression}% 
