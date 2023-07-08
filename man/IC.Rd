\name{IC}
\alias{IC}
\alias{AIC.gamlss}
\alias{GAIC}
\alias{extractAIC.gamlss}
\alias{GAIC.table}
\alias{GAIC.scaled}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{Gives the GAIC for a GAMLSS Object}

\description{
 The function \code{GAIC()} calculates the Generalised Akaike information criterion (GAIC) for a given penalty \code{k} for a fitted GAMLSS object. 
 
 The function \code{AIC.gamlss()} is the method associated with a GAMLSS object of the generic function \code{AIC()}. Note that \code{GAIC()} is a synonymous of the function \code{AIC.gamlss}.

 The function \code{IC()} is an old version of   \code{GAIC()}.

 The function \code{GAIC.table()} produces a table with different models as rows and different penalties, \code{k}, as columns.
 
The function \code{GAIC.scaled()} produces, [for a given set of different fitted models or for a table produced by \code{chooseDist()}], the  scaled  Akaike values (see Burnham  and  Anderson (2002) section 2.9 for a similar concept the GAIC weights. The scaled Akaike  should not be interpreted as posterior probabilities of models given the data but for model selection purpose they produce a scaled ranking of the model using their relative importance i.e. from the worst to the best model.   


 The function \code{extractAIC} is a the method associated a GAMLSS object of the generic function \code{extractAIC} and it is 
 mainly used in the \code{stepAIC} function. 
 
The function \code{Rsq} compute a generalisation of the R-squared for not normal models.
}
\usage{
IC(object, k = 2)
\method{AIC}{gamlss}(object, ..., k = 2, c = FALSE)
GAIC(object, ..., k = 2, c = FALSE )
GAIC.table(object, ..., k = c(2, 3.84, round(log(length(object$y)), 2)),
           text.to.show=NULL)
GAIC.scaled(object,..., k = 2, c = FALSE, plot = TRUE,   
              text.cex = 0.7, which = 1, diff.dev = 1000,
              text.to.show = NULL, col = NULL, horiz = FALSE) 
\method{extractAIC}{gamlss}(fit, scale, k = 2, c = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{an gamlss fitted model(s) [or for \code{GAIC.scaled()} a table 
  produced by \code{chooseDist()}].}
  \item{fit}{an gamlss fitted model}
  \item{\dots}{allows several GAMLSS object to be compared using a GAIC}
  \item{k}{the penalty with default \code{k=2.5}}
  \item{c}{whether the corrected AIC, i.e. AICc, should be used, note that it applies only when  \code{k=2}}
  \item{scale}{this argument is not used in gamlss}
    \item{plot}{whether to plot the ranking in \code{GAIC.scaled()}.}
  \item{text.cex}{the size of the models/families in the text of the plot of \code{GAIC.scaled()}.}
  \item{diff.dev}{this argument prevents models with a difference in deviance greater than \code{diff.dev} from the `best' model to be considered (or plotted).}
  \item{which}{which column of GAIC scaled to plot in \code{GAIC.scaled()}.}
  \item{text.to.show}{if NULL, \code{GAIC.scaled()} shows the model names otherwise the character in this list}
   \item{col}{The colour of the bars in  \code{GAIC.scaled()}}
   \item{horiz}{whether to plot the bars vertically (default) or horizontally}
   
}

\value{
  The function \code{IC()} returns the GAIC for given penalty k of the GAMLSS object.
  The function \code{AIC()} returns a matrix contains the df's and the GAIC's for given penalty k. 
  The function \code{GAIC()} returns identical results to \code{AIC}. 
  The function \code{GAIC.table()} returns a table which its rows showing different models and its columns different \code{k}'s.  
  The function \code{extractAIC()} returns vector of length two with the degrees of freedom and the AIC criterion.  
}
\references{ 

Burnham K. P.  and  Anderson D. R (2002). \emph{Model Selection and Multi model Inference
A Practical Information-Theoretic Approach}, Second Edition, Springer-Verlag New York, Inc.

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
\author{Mikis Stasinopoulos \email{d.stasinopoulos@londonmet.ac.uk} }

\seealso{ \code{\link{gamlss}} }
\examples{
data(abdom)
m1 <-  gamlss(y~x, family=NO, data=abdom)
IC(m1)
extractAIC(m1,k=2)
m2 <-  gamlss(y~x, sigma.fo=~x, family=NO, data=abdom)
m3 <-  gamlss(y~pb(x), sigma.fo=~x, family=NO, data=abdom)
m4 <-  gamlss(y~pb(x), sigma.fo=~pb(x), family=NO, data=abdom)
AIC(m1,m2, m3, m4)
AIC(m1,m2, m3, m4, c=TRUE)
AIC(m1,m2, m3, m4, k=3)
GAIC.table(m1,m2, m3, m4)
GAIC.scaled(m1,m2, m3, m4)
\dontrun{
MT <- chooseDist(m3)
GAIC.scaled(MT)
GAIC.scaled(MT, which=2)}
}

\keyword{regression}% 
