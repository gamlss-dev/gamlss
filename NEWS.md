# Version 5.4-19

* The package is now hosted on GitHub at
  <https://github.com/mstasinopoulos/GAMLSS-original/>.

* Add a new `prodist()` method for extracting fitted (in-sample) or predicted
  (out-of-sample) probability distributions from gamlss models (contributed by
  [Achim Zeileis](https://www.zeileis.org/)). This enables the workflow from the
  [distributions3](https://CRAN.R-project.org/package=distributions3) package for all
  distributions provided by `gamlss.dist`. The idea is that the `distributions3`
  objects encapsulate all information needed to obtain moments (mean, variance,
  etc.), probabilities, quantiles, etc. with a unified interface. See the
  [useR! 2022 presentation](https://www.zeileis.org/news/user2022/) by
  Zeileis, Lang, and Hayes for an overview.
