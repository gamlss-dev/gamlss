# create 23-6-18
#--------------------------------------------------------------
# TO DO
# i)   what action we should have with factors?
# ii)  what happents if linear terms? the first derivative look funny
# iii) check for different values of the other variables what happening  
# iv)  binomial or count data?
# -------------------------------------------------------------
getQuantile <- function (obj = NULL, 
                        term = NULL, 
                    quantile = .9,
                        data = NULL, 
                    n.points = 100, 
                  #parameter = c("mu", "sigma", "nu", "tau"), 
                     #  type = c("response", "link"), 
                         how = c("median", "last"), 
                    fixed.at = list(), 
                        plot = FALSE) 
{
  if (is.null(obj) || !class(obj)[1] == "gamlss") 
    stop("Supply a standard GAMLSS model in obj")
  if (is.null(term)) stop("The model term is not set")
   how <- match.arg(how)
  #type <- match.arg(type)
  #parameter <- match.arg(parameter)
  if (any(grepl("data", names(obj$call)))) {
    DaTa <- if (startsWith(as.character(obj$call["data"]), "na.omit")) 
      eval(parse(text = as.character(obj$call["data"])))
    else get(as.character(obj$call["data"]))
  }
  else if (is.null(data)) 
    stop("The data argument is needed in obj")
  mat <- matrix(0, nrow = dim(DaTa)[1] + n.points, ncol = dim(DaTa)[2])
         dat.temp <- as.data.frame(mat)
  names(dat.temp) <- v.names <- names(DaTa)
  pos <- which(names(dat.temp) == term)
  if (pos < 1) 
    stop("supply a continuous term")
# Think aboutfactors
  if (is.factor(DaTa[, pos])) 
    stop("the getQuantile() is not suitable for factors, what shall we Do??")
  xvar <- seq(from = min(DaTa[, pos]), to = max(DaTa[, pos]), 
              length.out = n.points)
  for (i in 1:dim(dat.temp)[2]) {
    if (pos == i) {
      dat.temp[, i] <- c(DaTa[, i], xvar)
    }
    else {
      ma <- fixed.at[[v.names[i]]]
      if (is.null(ma)) {
        if (how == "median") {
          ma <- if (is.factor(DaTa[, i])) 
            levels(DaTa[, i])[which.max(table(DaTa[, i]))]
          else median(DaTa[, i])
        }
        if (how == "last") {
          ma <- if (is.factor(DaTa[, i])) 
            levels(DaTa[, i])[which.max(table(DaTa[, 
                                                   i]))]
          else tail(DaTa[, i], 1)
        }
      }
      dat.temp[, i] <- c(DaTa[, i], rep(ma, n.points))
    }
  }
      pdf <- obj$family[1]
   binom  <- pdf%in%gamlss::.gamlss.bi.list # whether binomial
     qfun <- paste("q", obj$family[[1]],sep="")
     lpar <- eval(parse(text=pdf))()$nopar
  if (binom) {bd <- obj$bd ; Y <- obj$y}
    # get the prediction
       pp <-  predictAll(obj, newdata = tail(dat.temp, n.points), output="matrix")
  # get the quantile 
  # if (binom)
  #      {
  #      pp <-  predictAll(obj, newdata = tail(dat.temp, n.points), output="matrix")
  #        DevIncr <- switch(lpar, 
  #                          fn( Y[i], mu = pp[,"mu"], bd=bd[i]),   # 1
  #                          fn( Y[i], mu = pp[,"mu"],              # 2
  #                              sigma = pp[,"sigma"], bd=bd[i]),                        
  #                          fn( Y[i], mu = pp[,"mu"],              # 3
  #                              sigma = pp[,"sigma"],
  #                              nu = pp[,"nu"], bd=bd[i]),
  #                          fn( Y[i], mu = pp[,"mu"],              # 4
  #                              sigma = pp[,"sigma"],  
  #                              nu = pp[,"nu"],
  #                              tau = pp[,"tau"],bd=bd[i]))
  #      } else
  #      {
qq  <- switch(lpar, 
              eval(call(qfun, p= quantile, mu=pp[,"mu"])),       # 1
              eval(call(qfun, p= quantile, mu=pp[,"mu"], sigma=pp[,"sigma"])),        # 2
              eval(call(qfun, p= quantile, mu=pp[,"mu"], sigma=pp[,"sigma"],  nu=pp[,"nu"])),  # 3                   
              eval(call(qfun, p= quantile, mu=pp[,"mu"], sigma=pp[,"sigma"],  nu=pp[,"nu"], tau=pp[,"tau"])))
#       }  
     #  qq <- call(qfun, p= quantile, mu=pp[,"mu"], sigma=pp[,"sigma"])
                   
 theFun <- splinefun(xvar, qq)
  if (plot) {
    layout(matrix(c(1, 1, 2, 2), 2, 2, byrow = TRUE))
    Ylab <- paste("quantile=",as.character(quantile),sep="")
    plot(theFun(xvar) ~ xvar, ylab = Ylab, xlab = term, 
         type = "l")
    plot(theFun(xvar, deriv = 1) ~ xvar, xlab = term, ylab = "d/dx", 
         type = "l")
    abline(h = 0)
    layout(1)
  }
  invisible(theFun)
}