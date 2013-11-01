### function to plot additive terms in the fitting of a GAMLSS model
### it is based on the termplot() function in R
### but there were are several bugs corrected
### 1. if what parameters is constant an error is given
### 2. if interaction are present then a warning is given 
###    rather that crashing as is termplot() 
### Author: Mikis Stasinopoulos
### bugs: if the additve terms depends on more that one variable 
###       i.e. loess(x1,x2) produces rubish
term.plot <- function (object, 
                        what = c("mu","sigma","nu","tau"),  
                        data = NULL, 
                       envir = environment(formula(object)), 
               partial.resid = FALSE, 
                         rug = FALSE, 
                       terms = NULL, 
                          se = TRUE, 
                       xlabs = NULL, 
                       ylabs = NULL, 
                        main = NULL, 
                    col.term = 2, 
                    lwd.term = 1.5, 
                      col.se = "orange", 
                      lty.se = 2, 
                      lwd.se = 1, 
                     col.res = "gray",  
                         cex = 1, 
                         pch = par("pch"), 
                    col.smth = "darkred", 
                    lty.smth = 2, 
                   span.smth = 2/3, 
                         ask = interactive() && nb.fig < n.tms &&.Device != "postscript", 
           use.factor.levels = TRUE, 
                      smooth = NULL,  
                        ylim = "common", # 
                             ...) 
{
## only for gamlss objects 
    if (!is.gamlss(object))  stop(paste("This is not an gamlss object", "\n", "")) 
           what <- match.arg(what) 
## checking the parameter
    if (!what%in%object$par) stop(paste(what,"is not a parameter in the object","\n"))
## getting the specific term for plotting
    which.terms <- terms
## get all terms
      par.terms <- object[[paste(what, "terms", sep=".")]]
      par.attr  <- attributes(par.terms)
          terms <- if (is.null(terms)) lpred(object, what = what, type = "terms", se.fit = se)
                  else lpred(object, what = what, type = "terms", se.fit = se, terms = terms)
          n.tms <- ncol(tms <- as.matrix(if (se) terms$fit  else terms))
## if the parameters has only a constant fitted stop
     if (n.tms == 0)
          stop("The model for ",what, " has only the constant fitted") 
## data frame      
             mf <- model.frame(object, what = what) 
## I am not sure that we need this          
    if (is.null(data)) 
           data <- eval(object$call$data, envir)
    if (is.null(data)) 
           data <- mf
    if (NROW(tms) < NROW(data)) 
     {
       use.rows <- match(rownames(tms), rownames(data))
     }
    else use.rows <- NULL
## whether there are interactions in the model
            nmt <- colnames(tms)
   Interactions <- par.attr$order > 1 
## if interaction nmt has to change 
      if (any(Interactions))
       {
         nmt <- nmt[!Interactions]
         if (!se) 
          { terms <- terms[,nmt, drop = FALSE]
          }
         else 
          { 
              terms$fit <- terms$fit[,nmt,  drop = FALSE]  
           terms$se.fit <- terms$se.fit[,nmt,  drop = FALSE] 
          }
           n.tms <- ncol(tms <- as.matrix(if (se) terms$fit  else terms))
          # I am assuming that 'terms' will be used wisely here 
          warning("interactions have been taken out from the plots")
       }         
             cn <- parse(text = nmt)
    if (!is.null(smooth)) 
         smooth <- match.fun(smooth)
    if (is.null(ylabs)) 
          ylabs <- paste("Partial for", nmt)
    if (is.null(main)) 
           main <- ""
    else if (is.logical(main)) 
           main <- if (main)  deparse(object$call, 500)
                  else ""
                 else if (!is.character(main)) 
                   stop("`main' must be TRUE, FALSE, NULL or character (vector).")
          main <- rep(main, length = n.tms)
            pf <- envir
       carrier <- function(term) {
                                  if (length(term) > 1)   carrier(term[[2]])
                                  else eval(term, data, enclos = pf)
                                 }
  carrier.name <- function(term) {
                                  if (length(term) > 1) carrier.name(term[[2]])
                                  else as.character(term)
                                 }
    if (is.null(xlabs))
       { xlabs <- unlist(lapply(cn, carrier.name))
       }
    if (partial.resid || !is.null(smooth)) 
       {
          pres <- residuals(object, what = what, type="partial")
        if (!is.null(which.terms)) 
          pres <- pres[, which.terms, drop = FALSE]
        if (any(Interactions))
          pres <- pres[, nmt, drop = FALSE]
       }
        is.fac <- sapply(nmt, function(i) is.factor(mf[, i]))
      se.lines <- function(x, iy, i, ff = 2) 
                  {
                     tt <- ff * terms$se.fit[iy, i]
                    lines(x, tms[iy, i] + tt, lty = lty.se, lwd = lwd.se, 
                          col = col.se)
                    lines(x, tms[iy, i] - tt, lty = lty.se, lwd = lwd.se, 
                         col = col.se)
                  }
        nb.fig <- prod(par("mfcol"))
    if (ask) {
            op <- par(ask = TRUE)
             on.exit(par(op))
              }
    ylims <- ylim
    if (identical(ylims, "common")) 
    {
      ylims <- if (!se) 
        range(tms, na.rm = TRUE)
      else range(tms + 1.05 * 2 * terms$se.fit, tms - 1.05 * 
                   2 * terms$se.fit, na.rm = TRUE)
      if (partial.resid) 
        ylims <- range(ylims, pres, na.rm = TRUE)
      if (rug) 
        ylims[1L] <- ylims[1L] - 0.07 * diff(ylims)
    }
    for (i in 1:n.tms) 
    {
      # if we need differnt y limit for each variable
        if (identical(ylim, "free"))
          { 
               ylims <- range(tms[, i], na.rm = TRUE)
            if (se) 
               ylims <- range(ylims, tms[, i] + 1.05 * 2 * terms$se.fit[, 
                            i], tms[, i] - 1.05 * 2 * terms$se.fit[, i], 
                            na.rm = TRUE)
            if (partial.resid)  
              ylims <- range(ylims, pres[, i], na.rm = TRUE)
            if (rug) 
           ylims[1] <- ylims[1] - 0.07 * diff(ylims)
           }
        # if is a factor
        if (is.fac[i]) 
        {
            ff <- mf[, nmt[i]]
            if (!is.null(object$na.action)) 
                ff <- naresid(object$na.action, ff)
            ll <- levels(ff)
            xlims <- range(seq(along = ll)) + c(-0.5, 0.5)
            xx <- as.numeric(ff)
            if (rug) 
            {
                xlims[1] <- xlims[1] - 0.07 * diff(xlims)
                xlims[2] <- xlims[2] + 0.03 * diff(xlims)
            }
            plot(1, 0, type = "n", xlab = xlabs[i], ylab = ylabs[i], 
                xlim = xlims, ylim = ylims, main = main[i], xaxt = "n", 
                ...)
            if (use.factor.levels) 
                axis(1, at = seq(along = ll), labels = ll, ...)
            else axis(1)
            for (j in seq(along = ll)) 
            {
                ww <- which(ff == ll[j])[c(1, 1)]
                jf <- j + c(-0.4, 0.4)
                lines(jf, tms[ww, i], col = col.term, lwd = lwd.term, 
                  ...)
                if (se) 
                  se.lines(jf, iy = ww, i = i)
            }
        }  
        else 
        { # here is where changes had to be made at the moment every pass
          # cn is expression
            xx <- carrier(cn[[i]]) # ds Friday, October 9, 2009 at 13:13
          # why we need this ??? is it for random??
            if (is.factor(xx)) xx <- seq(along = levels(xx))
            if (!is.null(use.rows)) # in case soem of the rows are not used 
                xx <- xx[use.rows]
            xlims <- range(xx, na.rm = TRUE)
            if (rug) 
                xlims[1] <- xlims[1] - 0.07 * diff(xlims)
            oo <- order(xx)
            plot(xx[oo], tms[oo, i], type = "l", xlab = xlabs[i], 
                ylab = ylabs[i], xlim = xlims, ylim = ylims, 
                main = main[i], col = col.term, lwd = lwd.term, 
                ...)
            if (se) 
                se.lines(xx[oo], iy = oo, i = i)
        }
        if (partial.resid) {
            if (!is.fac[i] && !is.null(smooth)) {
                smooth(xx, pres[, i], lty = lty.smth, cex = cex, 
                  pch = pch, col = col.res, col.smooth = col.smth, 
                  span = span.smth)
            }
            else points(xx, pres[, i], cex = cex, pch = pch, 
                col = col.res)
        }
        if (rug) {
            n <- length(xx)
            lines(rep.int(jitter(xx), rep.int(3, n)), rep.int(ylims[1] + 
                c(0, 0.05, NA) * diff(ylims), n))
            if (partial.resid) 
                lines(rep.int(xlims[1] + c(0, 0.05, NA) * diff(xlims), 
                  n), rep.int(pres[, i], rep.int(3, n)))
        }
    }
    invisible(n.tms)
}







    
  
