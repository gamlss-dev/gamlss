\name{gamlss} 
\alias{gamlss}
\alias{is.gamlss}
\alias{gamlssNews}

\title{Generalized Additive Models for Location Scale and Shape}
\description{
  Returns an object of class "gamlss", which is a generalized additive model for location scale and shape (GAMLSS).
  The function \code{gamlss()} is very similar to the \code{gam()} function in S-plus (now also in R in package \code{gam}), but
  can fit more distributions (not only the ones belonging to the  exponential family) and can model all the parameters of the
 distribution as functions of the explanatory variables (e.g. using linear, non-linear, smoothing, loess and random effects terms).

This implementation of \code{gamlss()} allows modelling of up to four
parameters in a distribution family, which are conventionally called \code{mu}, \code{sigma}, \code{nu} and \code{tau}.

The function \code{gamlssNews()} shows what is new in the current implementation.
 }

\usage{
gamlss(formula = formula(data), sigma.formula = ~1, 
        nu.formula = ~1, tau.formula = ~1, family = NO(), 
        data = sys.parent(), weights = NULL, 
        contrasts = NULL, method = RS(),  start.from = NULL,  
        mu.start = NULL,  sigma.start = NULL, 
        nu.start = NULL, tau.start = NULL, 
        mu.fix = FALSE, sigma.fix = FALSE, nu.fix = FALSE, 
        tau.fix = FALSE, control = gamlss.control(...), 
        i.control = glim.control(...), ...)
is.gamlss(x)
gamlssNews()
}


\arguments{
  \item{formula}{a formula object, with the response on the left of an ~ operator, and the terms, separated by \eqn{+}{+} operators, on the right.
  Nonparametric smoothing terms are indicated by \code{pb()} for penalised beta splines, \code{cs} for smoothing splines, \code{lo} for loess smooth terms and \code{random} or \code{ra}
  for random terms, e.g. \code{y~cs(x,df=5)+x1+x2*x3}.  Additional smoothers can be added by creating the appropriate interface.
  Interactions with nonparametric smooth terms are not fully supported, but will not produce errors;
  they will simply produce the usual parametric interaction 
  }
  \item{sigma.formula}{a formula object for fitting a model to the sigma parameter, as in the formula above,
         e.g.  \code{sigma.formula=~cs(x,df=5)}. It can be abbreviated to \code{sigma.fo=~cs(x,df=5)}. }
  \item{nu.formula}{a formula object for fitting a model to the nu parameter, e.g. \code{nu.fo=~x} }
  \item{tau.formula}{a formula object for fitting a model to the tau parameter, e.g. \code{tau.fo=~cs(x,df=2)}}
  \item{family}{a \code{\link[gamlss.dist]{gamlss.family}} object, which is used to define the distribution and the link functions of the various parameters. 
   The distribution families supported by \code{gamlss()} can be found in \code{\link[gamlss.dist]{gamlss.family}}.
   Functions such as BI() (binomial) produce a family object. Also can be given without the parentheses i.e. \code{BI}.
   Family functions can take arguments, as in \code{BI(mu.link=probit)}} 
  \item{data}{a data frame containing the variables occurring in the formula. If this is missing, the variables should be on the search list.
     e.g. \code{data=aids}}
  %\item{na.action}{a function which indicates what should happen when the data contain \code{NA}s.  The default is set by the \code{na.action} setting
  %        of \code{options}, and is \code{na.omit}. Other options are \code{na.fail} and \code{na.exclude}. \code{na.omit} omits all the 
  %        rows in the data frame with the \code{NA}s.  \code{na.fail} will give an error if there are \code{NA}s in the data. 
  %        \code{na.exclude} behaves like \code{na.omit} as far the fitting of the model is concern but it gives correct length fitted values, 
  %        linear predictors and residuals with \code{NA}s in the missing rows}
  \item{weights}{a vector of weights. Note that this is not the same as in the glm() or gam() function. 
  Here weights can be used to weight out observations (like in \code{subset}) or for a weighted likelihood analysis where 
  the contribution of the observations to the likelihood differs according to \code{weights}. 
  The length of \code{weights} must be the same as the number of observations in the data. 
  By default, the weight is set to one. To set  weights to vector w use \code{weights=w}}
  \item{contrasts}{list of contrasts to be used for some or all of the factors appearing as variables in the model formula. 
  The names of the list should be the names of the corresponding variables. 
  The elements should either be contrast-type matrices (matrices with as many rows as levels of the factor and with columns 
  linearly independent of each other and of a column of ones), or else they should be functions that compute such contrast matrices.}
  \item{method}{the current algorithms for GAMLSS are RS(), CG() and mixed(). i.e. \code{method=RS()} will use the Rigby and Stasinopoulos algorithm, 
                 \code{method=CG()} will use the Cole and Green algorithm and 
  \code{mixed(2,10)} will use the RS algorithm twice before switching to the Cole and Green algorithm for up to 10 extra iterations}
  \item{start.from}{a fitted GAMLSS model which the fitted values will be used as staring values for the current model}
  \item{mu.start}{ vector or scalar of initial values for the location parameter mu e.g.  \code{mu.start=4}  }
  \item{sigma.start}{ vector or scalar of initial values for the scale parameter sigma e.g. \code{sigma.start=1} }
  \item{nu.start}{ vector or scalar of initial values for the parameter nu e.g. \code{nu.start=3} }
  \item{tau.start}{ vector or scalar of initial values for the location parameter tau e.g. \code{tau.start=2}}
  \item{mu.fix}{ whether the mu parameter should be kept fixed in the fitting processes e.g. \code{mu.fix=FALSE} }
  \item{sigma.fix}{whether the sigma parameter should be kept fixed in the fitting processes e.g. \code{sigma.fix=FALSE} }
  \item{nu.fix}{ whether the nu parameter should be kept fixed in the fitting processes e.g. \code{nu.fix=FALSE} }
  \item{tau.fix}{whether the tau parameter should be kept fixed in the fitting processes e.g. \code{tau.fix=FALSE}}
  \item{control}{this sets the control parameters of the outer iterations algorithm. The default setting is the \code{\link{gamlss.control}} function }
  \item{i.control}{this sets the control parameters of the inner iterations of the RS algorithm. The default setting is the \code{\link{glim.control}} function  }
  \item{\dots}{for extra arguments}
  \item{x}{an object}
} 

\details{The Generalized Additive Model for Location, Scale and Shape
 is a general class of statistical models for a univariate
response variable. The model assumes independent observations of the response variable
y given the parameters, the explanatory variables and the values
of the random effects. The distribution for the response variable
in the GAMLSS can be selected from a very general family of
distributions including highly skew and/or kurtotic continuous and
discrete distributions, see \code{\link[gamlss.dist]{gamlss.family}}. The systematic part of the model is
expanded to allow modelling not only of the mean (or location) parameter,
but also of the other parameters of the distribution of y, as
linear parametric and/or additive nonparametric  (smooth)
functions of explanatory variables and/or random effects terms.
Maximum (penalized) likelihood estimation is used to fit the
(non)parametric models. A Newton-Raphson/Fisher scoring algorithm
is used to maximize the (penalized) likelihood. The additive terms
in the model are fitted using a backfitting algorithm.


\code{is.gamlss} is a short version is \code{is(object,"gamlss")}

} 

\value{
  Returns a gamlss object with components
  \item{family}{the distribution family of the gamlss object (see \link[gamlss.dist]{gamlss.family})}
  \item{parameters }{the name of the fitted parameters i.e. \bold{\code{mu}}, \bold{\code{sigma}}, \bold{\code{nu}}, \bold{\code{tau}} }
  \item{call}{the call of the gamlss function}
  \item{y}{the response variable}
  \item{control}{the gamlss fit control settings}
  \item{weights}{the vector of weights}
  \item{G.deviance}{the global deviance}
  \item{N}{the number of observations in the fit}
  \item{rqres}{a function to calculate the normalized (randomized) quantile residuals of the object}
  \item{iter}{the number of external iterations in the fitting process}
  \item{type}{the type of the distribution or the response variable (continuous or discrete)}
  \item{method}{which algorithm is used for the fit, RS(), CG() or mixed()}
  \item{converged}{whether the model fitting has have converged}
  \item{residuals}{the normalized (randomized) quantile residuals of the model}
 \item{mu.fv}{the fitted values of the mu model, also  sigma.fv, nu.fv, tau.fv for the other parameters if present}
 \item{mu.lp}{the linear predictor of the mu model, also  sigma.lp, nu.lp, tau.lp for the other parameters if present}
 \item{mu.wv}{the working variable of the mu model, also  sigma.wv, nu.wv, tau.wv for the other parameters if present}
 \item{mu.wt}{the working weights of the mu model, also  sigma.wt, nu.wt, tau.wt for the other parameters if present}
 \item{mu.link}{the link function for the mu model, also  sigma.link, nu.link, tau.link for the other parameters if present}
 \item{mu.terms}{the terms for the mu model, also  sigma.terms, nu.terms, tau.terms for the other parameters if present}
 \item{mu.x}{the design matrix for the mu, also  sigma.x, nu.x, tau.x for the other parameters if present}
 \item{mu.qr}{the QR decomposition of the mu model, also sigma.qr, nu.qr, tau.qr for the other parameters if present}
 \item{mu.coefficients}{the linear coefficients of the mu model, also  sigma.coefficients, nu.coefficients, tau.coefficients for the other parameters if present }
 \item{mu.formula}{the formula for the mu model, also  sigma.formula, nu.formula, tau.formula for the other parameters if present}
 \item{mu.df}{the mu degrees of freedom also  sigma.df, nu.df, tau.df for the other parameters if present}
 \item{mu.nl.df}{the non linear degrees of freedom, also sigma.nl.df, nu.nl.df, tau.nl.df for the other parameters if present} 
 \item{df.fit}{the total degrees of freedom use by the model}
 \item{df.residual}{the residual degrees of freedom left after the model is fitted}
 \item{aic}{the Akaike information criterion }
 \item{sbc}{the Bayesian information criterion }
 
} 

\references{

Rigby, R. A. and  Stasinopoulos D. M. (2005). Generalized additive models for location, scale and shape,(with discussion), 
\emph{Appl. Statist.}, \bold{54}, part 3, pp 507-554.

Stasinopoulos D. M., Rigby R.A. and Akantziliotou C. (2006) Instructions on how to use the GAMLSS package in R.
Accompanying documentation in the current GAMLSS  help files, (see also  \url{http://www.gamlss.org/}).

Stasinopoulos D. M. Rigby R.A. (2007) Generalized additive models for location scale and shape (GAMLSS) in R.
\emph{Journal of Statistical Software}, Vol. \bold{23}, Issue 7, Dec 2007, \url{http://www.jstatsoft.org/v23/i07}.
}


\author{ Mikis Stasinopoulos \email{d.stasinopoulos@londonmet.ac.uk}, Bob Rigby \email{r.rigby@londonmet.ac.uk}, Calliope Akantziliotou and Vlasios Voudouris \email{V.Voudouris@londonmet.ac.uk}}

\note{The following generic functions can be used with a GAMLSS object: \code{print}, \code{summary},  \code{fitted},  \code{coef},  
          \code{residuals}, \code{update},  \code{plot}, \code{deviance}, \code{formula}   }

\section{Warning }{Respect the parameter hierarchy when you are fitting a model. For example a good model for \code{mu} should be fitted before a model for \code{sigma} is fitted} 

\seealso{ \code{\link[gamlss.dist]{gamlss.family}}, \code{\link{pdf.plot}}, \code{\link{find.hyper}}  }

\examples{ 
data(abdom)
mod<-gamlss(y~pb(x),sigma.fo=~pb(x),family=BCT, data=abdom, method=mixed(1,20))
plot(mod)
rm(mod)
}


\keyword{regression}% 
