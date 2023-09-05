\name{prodist.gamlss}

\alias{prodist.gamlss}

\title{Extracting Fitted or Predicted Probability Distributions from gamlss Models}

\description{
Methods for \pkg{gamlss} model objects for extracting fitted (in-sample) or
predicted (out-of-sample) probability distributions as \pkg{distributions3}
objects.
}

\usage{
\method{prodist}{gamlss}(object, ...)
}
\arguments{
  \item{object}{A model object of class \code{\link{gamlss}}.}
  \item{...}{Arguments passed on to \code{\link{predictAll}}, e.g., \code{newdata}.}
}

\details{
To facilitate making probabilistic forecasts based on \code{\link{gamlss}}
model objects, the \code{\link[distributions3]{prodist}} method extracts fitted
or predicted probability \code{distribution} objects. Internally, the
\code{\link{predictAll}} method is used first to obtain the distribution
parameters (\code{mu}, \code{sigma}, \code{tau}, \code{nu}, or a subset thereof).
Subsequently, the corresponding \code{distribution} object is set up using the
\code{\link[gamlss.dist]{GAMLSS}} class from the \pkg{gamlss.dist} package,
enabling the workflow provided by the \pkg{distributions3} package (see Zeileis
et al. 2022).

Note that these probability distributions only reflect the random variation in
the dependent variable based on the model employed (and its associated
distributional assumption for the dependent variable). This does not capture the
uncertainty in the parameter estimates.
}

\value{
An object of class \code{GAMLSS} inheriting from \code{distribution}.
}

\references{
Zeileis A, Lang MN, Hayes A (2022).
\dQuote{distributions3: From Basic Probability to Probabilistic Regression.}
Presented at \emph{useR! 2022 - The R User Conference}.
Slides, video, vignette, code at \url{https://www.zeileis.org/news/user2022/}.
}

\seealso{
\code{\link[gamlss.dist]{GAMLSS}}, \code{\link{predictAll}}
}

\examples{
\dontshow{ if(!requireNamespace("distributions3")) {
  if(interactive() || is.na(Sys.getenv("_R_CHECK_PACKAGE_NAME_", NA))) {
    stop("not all packages required for the example are installed")
  } else q() }
}
## packages, code, and data
library("gamlss")
library("distributions3")
data("cars", package = "datasets")

## fit heteroscedastic normal GAMLSS model
## stopping distance (ft) explained by speed (mph)
m <- gamlss(dist ~ pb(speed), ~ pb(speed), data = cars, family = "NO")

## obtain predicted distributions for three levels of speed
d <- prodist(m, newdata = data.frame(speed = c(10, 20, 30)))
print(d)

## obtain quantiles (works the same for any distribution object 'd' !)
quantile(d, 0.5)
quantile(d, c(0.05, 0.5, 0.95), elementwise = FALSE)
quantile(d, c(0.05, 0.5, 0.95), elementwise = TRUE)

## visualization
plot(dist ~ speed, data = cars)
nd <- data.frame(speed = 0:240/4)
nd$dist <- prodist(m, newdata = nd)
nd$fit <- quantile(nd$dist, c(0.05, 0.5, 0.95))
matplot(nd$speed, nd$fit, type = "l", lty = 1, col = "slategray", add = TRUE)

## moments
mean(d)
variance(d)

## simulate random numbers
random(d, 5)

## density and distribution
pdf(d, 50 * -2:2)
cdf(d, 50 * -2:2)

## Poisson example
data("FIFA2018", package = "distributions3")
m2 <- gamlss(goals ~ pb(difference), data = FIFA2018, family = "PO")
d2 <- prodist(m2, newdata = data.frame(difference = 0))
print(d2)
quantile(d2, c(0.05, 0.5, 0.95))

## note that log_pdf() can replicate logLik() value
sum(log_pdf(prodist(m2), FIFA2018$goals))
logLik(m2)
}

\keyword{distribution}
