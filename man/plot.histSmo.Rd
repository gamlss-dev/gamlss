\name{plot.histSmo}
\alias{plot.histSmo}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
A Plotting Function for density estimator object histSmo 
}
\description{
Plots the estimated density or its c.d.f function or its inverse cdf function  
}
\usage{
\method{plot}{histSmo}(x, type = c("hist", "cdf", "invcdf"), ...)
}

\arguments{
  \item{x}{An histSmo object}
  \item{type}{Different plots: a histogram and density estimator, a cdf function or an inverse cdf function. 
}
  \item{\dots}{
for further arguments
}
}

\value{
returns the relevant plot
}
\references{
Eilers, P. (2003). A perfect smoother. \emph{Analytical Chemistry}, 75: 3631-3636.

Eilers, P. H. C. and Marx, B. D. (1996). Flexible smoothing with
B-splines and penalties (with comments and rejoinder). \emph{Statist. Sci},
\bold{11}, 89-121.

Lindsey, J.K. (1997) \emph{Applying Generalized Linear Models}. New York: Springer-Verlag. 
ISBN 0-387-98218-3 


Rigby, R. A. and  Stasinopoulos D. M. (2005). Generalized additive models for location, scale and shape,(with discussion), 
\emph{Appl. Statist.}, \bold{54}, part 3, pp 507-554.

Stasinopoulos D. M. Rigby R.A. (2007) Generalized additive models for location scale and shape (GAMLSS) in R.
\emph{Journal of Statistical Software}, Vol. \bold{23}, Issue 7, Dec 2007, \url{http://www.jstatsoft.org/v23/i07}.
}
\author{
Mikis Stasinopoulos, Paul Eilers, Bob Rigby, Vlasios Voudouris and Majid  Djennad
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{histSmo}}
}
\examples{
Y <- rPARETO2(1000) 
m1<- histSmo(Y, lower=0, save=TRUE)
plot(m1)
plot(m1, "cdf")
plot(m1, "invcdf")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{distribution}
