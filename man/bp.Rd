\name{bp}
\alias{bp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Bucket plot}
\description{
A bucket plot is a graphical way to check the skewness and kurtosis of a continuous variable or the residuals of a fitted GAMLSS model. It plots the transformed moment skewness and  transformed moment kurtosis of the variable (or residuals) together with a cloud of points obtained using a non-parametric bootstrap from the original variable (or residuals). It also provides a graphical way of performing the  Jarque-Bera test (JarqueandBera,1980). 

There are two different bucket plots specified by the \code{type} argument:

i) the \code{moment} bucket   and 
ii) the \code{centile} bucket which itself can be  \code{central} or \code{tail} one.

}
\usage{
bp(obj = NULL, weights = NULL, 
      type = c("moment", "centile.central", "centile.tail"), 
      xvar = NULL, bootstrap = TRUE, no.bootstrap = 99, 
      col.bootstrap = c("lightblue", "pink", "khaki", 
      "thistle", "tan", "sienna1","steelblue", "coral", "gold", 
      "cyan"), 
      pch.bootstrap = rep(21, 10), asCharacter = TRUE, 
      col.point = rep("black", 10), pch.point = 1:10, 
      lwd.point = 2, text.to.show = NULL, cex.text = 1.5, 
      col.text = "black", show.legend = FALSE, n.inter = 4, 
      xcut.points = NULL, overlap = 0, show.given = TRUE, 
      cex = 1, pch = 21, data = NULL, 
      bar.bg = c(num = "lightblue", fac = "pink"), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obj}{A \code{gamlss} fitted object.}
  \item{weights}{prior weights.}
  \item{type}{type of bucket plot whether "moment", "centile.central", or "centile.tail".}
  \item{xvar}{the x-variable if need to split the bucket plot.}
  \item{bootstrap}{whether to bootstrap the skewness and kurtosis points}
  \item{no.bootstrap}{the number of the  bootstrap samples in the plot}
  \item{col.bootstrap}{the colour of the  bootstrap samples in the plot}
  \item{pch.bootstrap}{the character plotting symbol.}
  \item{asCharacter}{whether to plot the skewness and kurtosis as character or just  points.}
  \item{col.point}{the colout of the point is plotted as point}
  \item{pch.point}{the character symbol for the point}
  \item{lwd.point}{the width of the symbol}
  \item{text.to.show}{whether to show character for the model}
  \item{cex.text}{the \code{cex} of the text}
  \item{col.text}{the \code{colour} of the text}
  \item{show.legend}{whether to show the legend}
  \item{n.inter}{number of intervals}
  \item{xcut.points}{cut points for the \code{xvar} if need}
  \item{overlap}{whether the interval id \code{xvar} is set should overlap}
  \item{show.given}{showing the top part of the plot}
  \item{cex}{the \code{cex}}
  \item{pch}{the point character \code{pch}}
  \item{data}{if data has to be set}
  \item{bar.bg}{the backgroud color of the bars in the top of the figure}
  \item{\dots}{other arguments}
}

\value{
A plot displaying the transformed moment skewness and  transformed moment kurtosis of the sample or residual of a model.
}
\references{
De Bastiani, F.,   Stasinopoulos, D. M., Rigby, R. A., Heller, G. Z., and Lucas A. (2022) Bucket Plot: A Visual Tool for Skewness and Kurtosis Comparisons. send for publication.

Rigby, R. A. and  Stasinopoulos D. M. (2005). Generalized additive models for location, scale and shape,(with discussion), 
\emph{Appl. Statist.}, \bold{54}, part 3, pp 507-554.

Rigby, R. A., Stasinopoulos, D. M.,  Heller, G. Z.,  and De Bastiani, F. (2019)
	\emph{Distributions for modeling location, scale, and shape: Using GAMLSS in R}, Chapman and Hall/CRC. \doi{10.1201/9780429298547}
	An older version can be found in \url{https://www.gamlss.com/}.

Stasinopoulos D. M. Rigby R.A. (2007) Generalized additive models for location scale and shape (GAMLSS) in R.
\emph{Journal of Statistical Software}, Vol. \bold{23}, Issue 7, Dec 2007, \doi{10.18637/jss.v023.i07}.

Stasinopoulos D. M., Rigby R.A., Heller G., Voudouris V., and De Bastiani F., (2017)
\emph{Flexible Regression and Smoothing: Using GAMLSS in R},  Chapman and Hall/CRC. \doi{10.1201/b21973}

Stasinopoulos, M. D., Rigby, R. A., and De Bastiani F., (2018) GAMLSS: a distributional regression approach, \emph{Statistical Modelling}, Vol. \bold{18}, 
pp, 248-273, SAGE Publications Sage India: New Delhi, India. \doi{10.1177/1471082X18759144}

(see also  \url{https://www.gamlss.com/}).
}
\author{Mikis Stasinopoulos, \email{d.stasinopoulos@londonmet.ac.uk}, Bob Rigby 
            \email{r.rigby@londonmet.ac.uk} and Fernanda De Bastiani}
\note{
The  bucket plot provides an additional residual diagnostic tool that can be used for fitted model checking, alongside other diagnostic tools, for example   worm plots, and Q (and Z) statistics.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{\link[gamlss]{wp}, \link[gamlss]{Q.stats}}}
\examples{
m1 <- gamlss(R~pb(Fl)+pb(A), data=rent, family=GA)
bp(m1)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
\keyword{regression}
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
\concept{diagnostic}
% Use only one concept per line.
