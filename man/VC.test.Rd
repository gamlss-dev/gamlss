\name{VC.test}
\alias{VC.test}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Vuong and Clarke tests
}
\description{
The Vuong and Clarke tests for GAMLSS fitted models.
}
\usage{
VC.test(obj1, obj2, sig.lev = 0.05)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obj1}{The first  fitted gamlss object}
  \item{obj2}{The second  fitted gamlss object}
  \item{sig.lev}{Significance level used for testing.}
  }
\details{
The Vuong (1989) and Clarke (2007) tests are likelihood-ratio-based tests for model selection that use the Kullback-Leibler information criterion. 
The implemented tests can be used for choosing between two bivariate models which are non necessary nested.

In the Vuong test, the null hypothesis is that the two models are equally close to the actual model, whereas the alternative is that one model is closer. The test follows asymptotically a standard normal distribution under the null. Assume that the critical region is \code{(-c,c)}, where \code{c} is typically set to 1.96. If the value of the test is greater than \code{c} then we reject the null hypothesis that the models are equivalent in favour of the model in \code{obj1}. Vice-versa if the value is smaller than \code{-c} we reject the null hypothesis that the models are equivalent in favour of the model in \code{obj2}. If the value falls within \code{(-c,c0)} then we cannot discriminate between the two competing models given the data.

In the Clarke test, if the two models are statistically equivalent then the log-likelihood ratios of the observations should be evenly distributed around zero and around half of the ratios should be larger than zero. The test follows asymptotically a binomial distribution with parameters n and 0.5. Critical values can be obtained as shown in Clarke (2007). Intuitively, the model in \code{obj1} is preferred over that in \code{obj2} if the value of the test is significantly larger than its expected value under the null hypothesis ('code{n/2}), and vice versa. If the value is not significantly different from \code{n/2} then \code{obj1} can be thought of as equivalent to \code{obj2}.
}
\value{
For the Vuong test it returns its value and the decision and for the Clarke test returns the value the p-value and the decision. Decisions criteria are as discussed above.  
}
\references{
Clarke K. (2007), A Simple Distribution-Free Test for Non-Nested Model Selection. \emph{Political Analysis}, 15, 347-363.

Vuong Q.H. (1989), Likelihood Ratio Tests for Model Selection and Non-Nested Hypotheses. \emph{Econometrica}, 57(2), 307-333.

}
\author{
Mikis Stasinopoulos and Giampierro Marra
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{LR.test}}
}
\examples{
library(gamlss)
# fitting different models 
m0 <- gamlss(y~x+qrt, data=aids, family=PO)
m1 <- gamlss(y~pb(x)+qrt, data=aids, family=PO)
m2 <- gamlss(y~pb(x)+qrt, data=aids, family=NBI)
# comparison of the mdels
VC.test(m0,m2)
VC.test(m0,m1)
VC.test(m1,m2)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{regression}
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
