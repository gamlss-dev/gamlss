\name{plot2way}
\alias{plot2way}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to plot two interaction in a GAMLSS model
}
\description{
This function is designed to plot a factor to factor interaction in a GAMLSS model.
}
\usage{
plot2way(obj, terms = list(), what = c("mu", "sigma", "nu", "tau"), 
        show.legend = TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obj}{A gamlss model
}
  \item{terms}{this should be a character vector with the names of the two factors to be plotted}
  \item{what}{which parameters? \code{mu}, \code{sigma}, \code{nu}, or \code{tau} }
  \item{show.legend}{
whether to show the legend in the two way plot 
}
  \item{\dots}{Further arguments
}
}
\details{
This is an experimental function which should be use with prudence since no other check is done on whether this interaction interfere with other terms in the model
}
\value{
The function creates a 2 way interaction plot 
}
\references{
Rigby, R. A. and  Stasinopoulos D. M. (2005). Generalized additive models for location, scale and shape,(with discussion), 
\emph{Appl. Statist.}, \bold{54}, part 3, pp 507-554.

Stasinopoulos D. M., Rigby R.A. and Akantziliotou C. (2006) Instructions on how to use the GAMLSS package in R.
Accompanying documentation in the current GAMLSS  help files, (see also  \url{http://www.gamlss.org/}).  

Stasinopoulos D. M. Rigby R.A. (2007) Generalized additive models for location scale and shape (GAMLSS) in R.
\emph{Journal of Statistical Software}, Vol. \bold{23}, Issue 7, Dec 2007, \url{http://www.jstatsoft.org/v23/i07}.
}
\author{
Mikis Stasinopoulos
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{term.plot}},
}
\examples{
data(aids)
ti <- factor(c(rep(1,18),rep(2,27))) 
m1 <- gamlss(y~x+qrt*ti, data=aids, family=NBI)
m2 <- gamlss(y~x+qrt*ti, data=aids, family=NO)
plot2way(m1, c("qrt","ti")) 
plot2way(m1, c("ti", "qrt"))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{regression}
