# Latest mofifictions 
#  i) the average version of the rqres.plot() is modefified on the 21-Sep-2017
# ii) the function get.rqres() is introduced  
#-------------------------------------------------------------------
#-------------------------------------------------------------------
# new fuction for a lot of randomised 
get.rqres <- function(obj = NULL, howmany=10, order=FALSE) #, pit=FALSE
{  if (!is.gamlss(obj))  stop(paste("This is not an gamlss object", "\n", ""))
  # if (obj$type!="Discrete" ) stop(paste("This is not discrete distribution ", "\n", ""))
  # if (howmany>10&plot=="all")  stop(paste("You can only have 10 or less plots" , "\n", ""))
  # if (pit) obj$rqres[[1]][[length(obj$rqres)+1]] <- alist(pit=TRUE)
    w <- obj$weights
  if (all(trunc(w)==w)) # if frequenies as weights
  { 
    y  <- rep(obj$y, w)
    mu <- rep(fitted(obj, "mu"),w)
    if(any(obj$family%in%.gamlss.bi.list)){ bd <- rep(obj$bd,w)} # MS Wednesday, July 23, 2003 at 12:03   
    if ("sigma"%in%obj$parameters)  sigma <- rep(fitted(obj,"sigma"),w)
    if ("nu"%in%obj$parameters)        nu <- rep(fitted(obj,   "nu"),w)
    if ("tau"%in%obj$parameters)      tau <- rep(fitted(obj,  "tau"),w)  
  }
  else   # note that weights=1 and weights not frequencies are treated equal here and this could create problems in the future
  {
    y  <- obj$y
    mu <- fitted(obj)
    if(any(obj$family%in%.gamlss.bi.list)){ bd <- obj$bd} # MS Wednesday, July 23, 2003 at 12:03   
    if ("sigma"%in%obj$parameters)  sigma <- fitted(obj,"sigma")
    if ("nu"%in%obj$parameters)        nu <- fitted(obj,"nu")
    if ("tau"%in%obj$parameters)      tau <- fitted(obj,"tau")
  }
    rs <- matrix(0,ncol=howmany, nrow =length(y) ) 
   for (i in 1:howmany)
   {
         res <-  eval(obj$rqres)
     rs[,i]  <- if (order) res[order(res)] else res
   }
 rs
}
#-------------------------------------------------------------------
#-------------------------------------------------------------------
 rqres.plot <- function (obj = NULL, 
                     howmany = 6, 
                   plot.type = c("few", "all"),# 
                        type = c("wp", "QQ"),
                      #    cb = c(0.025, 0.975),
                        xlim = NULL, 
                        ylim = NULL,
                        ...)
{
  plot.type  <- match.arg(plot.type)
        type <- match.arg(type)
        var1 <- ceiling(howmany/2)
        lobj <- obj$noObs# length(fitted(obj))
         rs <- get.rqres(obj=obj, howmany=howmany, order=TRUE)
 switch(plot.type,
        "few"=  # only few plots
        {
        if (howmany >8) stop("howmany for \'few\' should be less than 8")  
        op <- par(mfrow=c(var1,2), col.axis="blue4", col.main="blue4", col.lab="blue4",col="darkgreen", bg="beige")
         on.exit(par(op))
           for (i in 1:howmany)
             {
                   res <-  rs[,i]
               if (type=="QQ"){#  ----  QQ plots here
                     rs1 <-    qqnorm(res, main = "Normal Q-Q Plot",
                                      xlab = "Theoretical Quantiles",
                                      ylab = "Sample Quantiles", 
                                      plot.it = TRUE, 
                                      col = "darkgreen")$y
                     lines(rs1, rs1, col="red" , lwd=.4, cex=.4 )     
                           } 
              if (type=="wp"){# ---- worm plots  plots here
                rs[,i] <- res
                wp(resid=res, ...)
              }     
            }
         },
        "all"=  # --- plotting all randomise quantile residuals 
        {
           rmedian <- apply(rs, MARGIN=1, "median")
          if (type=="QQ"){ # ----   QQ-plots
            x <-qqnorm(rs[,1],  plot=FALSE)$x
            plot(rs[,1]~x,    xlab = "Theoretical Quantiles",
                 ylab = "Sample Quantiles", type="n")
            for (i in 1:howmany)
            {
              points(rs[,i]~x, col="lightgray", pch=21)
            }
            points(rmedian~x,  pch=20, cex=.4 )
            lines(rmedian, rmedian, col="red" , lwd=.4, cex=.4 )
            }
          if (type=="wp"){# --- worm plots
                  x <- qqnorm(rmedian,  plot=FALSE)$x
                 yy <- rmedian - x
              if (is.null(xlim)) 
                   xlim <- range(x)
              if (is.null(ylim)) 
              {rr<-range(rs-x); ylim <- c(-max(rr),max(rr))} 
                
          plot(yy~x,ylab="Deviation", xlab="Unit normal quantile", 
               col="lightgray", type="n", ylim=ylim, xlim=xlim)
            for (i in 1:howmany)
            {
              yy <- rs[,i] - x
              points(yy~x, col="lightgray", pch=21)
            }
          yy <- rmedian - x
          points(yy~x,  pch=20, cex=.4 )
          grid(lty = "solid")
          abline(0, 0, lty = 2, col = gray(.6))
          abline(0, 100000, lty = 2, col = gray(.6))
          yuplim <- 10*sqrt(1/length(rmedian))
           level <- .95    
               z <- seq(-5,5,0.25)
               p <- pnorm(z)   
              se <- (1/dnorm(z))*(sqrt(p*(1-p)/length(rmedian)))     
             low <- qnorm((1-level)/2)*se
            high <- -low      
            lines(z, low, lty=2)
             lines(z, high, lty=2)
          }
    #  }, 
     # "cb"=
     # {
     #   rmedian <- apply(rs, MARGIN=1, "median")
     #     low <- apply(rs, MARGIN=1, "quantile",0.025)
     #   upper <- apply(rs, MARGIN=1, "quantile",0.975)
     #   if (type=="QQ"){
     #     x <-qqnorm(rs[,1],  plot=FALSE)$x
     #     plot(rmedian~x,    xlab = "Theoretical Quantiles",
     #          ylab = "Sample Quantiles", type="n")
     #     
     #     xx <- c(x,rev(x))
     #     yy <- c(low, upper)
     #     polygon(xx, yy, col = "lightgray", border = "darkgray")
     #     points(rmedian~x, pch=20, cex=.4 )
     #     lines(rmedian, rmedian, col="red" , lwd=.4, cex=.4 )
     #   }
     #   if (type=="wp"){
     #     # wp(resid=rmedian, ...)
     #     #----  
     #     x <- qqnorm(rmedian,  plot=FALSE)$x
     #     yy <- rmedian - x
     #     ylim.a <- max(yy)*2
     #     plot(yy~x,ylab="Deviation", xlab="Unit normal quantile", 
     #          col="lightgray", type="n", ylim=c(-ylim.a, ylim.a))
     #     for (i in 1:howmany)
     #     {
     #       yy <- rs[,i] - x
     #       points(yy~x, col="lightgray", pch=21)
     #     }
     #     yy <- rmedian - x
     #     points(yy~x, , pch=20, cex=.4 )
     #     grid(lty = "solid")
     #     abline(0, 0, lty = 2, col = 2)
     #     abline(0, 100000, lty = 2, col = 2)
     #     yuplim <- 10*sqrt(1/length(rmedian))
     #     level <- .95    
     #     z <- seq(-5,5,0.25)
     #     p <- pnorm(z)   
     #     se <- (1/dnorm(z))*(sqrt(p*(1-p)/length(rmedian)))     
     #     low <- qnorm((1-level)/2)*se
     #     high <- -low      
     #     lines(z, low, lty=2)
     #     lines(z, high, lty=2)
      # }
     })
     invisible(rs)      
}
