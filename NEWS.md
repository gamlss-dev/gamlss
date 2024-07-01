# Introduction

gamlss is a R package implementing the Generalised Additive Models for
 The package `gamlss` is an implementation of  Rigby and Stasinopoulos (2005),  Appl. Statist., 54,  pp. 507-554.

There are three book available for information; 

 1) "Flexible Regression and Smoothing: Using GAMLSS in R" 
explaining how the models can be used in R.

2) 
"Distributions for modeling location, scale and shape: Using GAMLSS in R" 
explaining the explicit and generated distributions available in the 
package gamlss.dist  

3)  
"Generized Additive Models for Location Scale and Shape: A distributional 
regression  approach with applications" 
explaining the different method for fitting GAMLSS i.e. penalised Likelihood, Bayesian and Boosting.  
 
More more information about books and papers related to GAMLSS can be found in
<https://www.gamlss.com/>.
 
 
The GitHub repository is now hosted under the new `gamlss-dev` organization:
  <https://github.com/gamlss-dev/gamlss/>.



# Version 5.4-23

- Tim Cole's suggestion in  `predictAll()` is added. This is to deal with the problem when `mu` is fixed.

- Tim Cole's suggestion in `summary()`is added. This to fix the problem when y~0, (that is, when there are no df's), to be incorporated in the `summary.gamlss()`.



# Version 5.4-21

* `predict()` do not  print the message "new prediction"

* `stepGAIC()` produce less lines in the output 




# Version 5.4-20

* The package is now hosted on GitHub at
  <https://github.com/gamlss-dev/gamlss/>.

* Add a new `prodist()` method for extracting fitted (in-sample) or predicted
  (out-of-sample) probability distributions from gamlss models (contributed by
  [Achim Zeileis](https://www.zeileis.org/)). This enables the workflow from the
  [distributions3](https://CRAN.R-project.org/package=distributions3) package for all
  distributions provided by `gamlss.dist`. The idea is that the `distributions3`
  objects encapsulate all information needed to obtain moments (mean, variance,
  etc.), probabilities, quantiles, etc. with a unified interface. See the
  [useR! 2022 presentation](https://www.zeileis.org/news/user2022/) by
  Zeileis, Lang, and Hayes for an overview.
