#------------------------------------------------------------------------------------------
#                                   pdf.plot
#                  Kalliope Akantziliotou, Mikis Stasinopoulos
#                         Thursday, February 6, 2003 at 12:04
#             last May 2019 Mikis and Fernanda
#------------------------------------------------------------------------------------------
#------------------------------------------------------
pdf.plot<- function (obj = NULL,
                     obs = c(1),
                  family = NO(),
                      mu = NULL,
                   sigma = NULL,
                      nu = NULL,
                     tau = NULL,
                    from = 0,
                      to = 10,
                     min = NULL,
                     max = NULL,
               no.points = 201, 
                no.title = FALSE,
                     col = gray(.4),
              y.axis.lim = 1.1, 
              frame.plot = TRUE,
                           ...)
{
  if (is.null(min)) {min <- from}
  if (is.null(max)) {max <- to}
gamlss.bi.list <- .binom
  if (!is.null(obj)) family <- obj[["family"]][1]
        fname <- if (is.name(family)) as.character(family)
  else if (is.character(family)) family
  else if (is.call(family)) as.character(family[[1]])
  else if (is.function(family)) deparse(substitute(family))
  else if (is(family, "gamlss.family"))  family$family[1]
  else stop("the family must be a character or a gamlss.family name")
    fam1 <- eval(parse(text=fname)) # the family to output
     fam <- as.gamlss.family(family) # this is created so I can get things
  dorfun <- paste("d",fname,sep="") # say dNO
   nopar <- fam$nopar # or fam1$nopar
    type <- fam$type
#--------------------------------------------------------- 
   y.var <- if(type=="Discrete")  seq(min,max,by=1)
            else seq(min,max,length=no.points)
  if(!is.null(obj))
  {
if (!is.gamlss(obj))  stop(paste("This is not an gamlss object", "\n", ""))
    lobs <- length(obs)
if (lobs > 9)   stop(paste("Use up to eight observations for ploting"))
  plots <- list(c(1,1),c(2,1),c(3,1),c(2,2),c(3,2),c(3,2),c(4,2),c(4,2))
    vvv <- unlist(plots[lobs])
if (type=="Discrete")   {typelh <- "h" } else   {typelh <- "l" }
  linter <- length(y.var)
if (any(fname%in%.gamlss.bi.list))
    {
      stop(paste("This function is not working for binomial type distibutions", "\n", ""))
    }
if ("mu"%in%names(fam$parameters))
    { 
  mu.var <- fitted(obj)[obs]
  if (!fam$mu.valid(mu.var))  stop( "`mu' parameter out of range")
    }
if ("sigma"%in%names(fam$parameters))
    {  
sigma.var <- fitted(obj,"sigma")[obs]
  if (!fam$sigma.valid(sigma.var))  stop( "`sigma' parameter out of range")
    }
if ("nu"%in%names(fam$parameters))
    {  
 nu.var <- fitted(obj,"nu")[obs]
  if (!fam$nu.valid(nu.var))  stop( "`nu' parameter out of range")
    }
if ("tau"%in%names(fam$parameters))
    { 
 tau.var <- fitted(obj,"tau")[obs]
 if (!fam$tau.valid(tau.var))  stop( "`tau' parameter out of range")
    }
     pdfArr <- array(0:0, c(linter,lobs))
title.label <- rep(NA, lobs) 
op <- par(mfrow=vvv, mar=par("mar")+c(0,0,0,0),
              col.axis=gray(.2),col.main=gray(.2),
              col.lab=gray(.2), cex=0.6, font=2)
# loop
for (j in 1:lobs)
 {
    pdf11 <-  switch (nopar,
                paste0("d",fname,"(y.var,mu.var[j])"),
                paste0("d",fname,"(y.var,mu.var[j],sigma.var[j])"),
                paste0("d",fname,"(y.var,mu.var[j],sigma.var[j], nu.var[j])"),
                paste0("d",fname,"(y.var,mu.var[j],sigma.var[j], nu.var[j], 
                       tau.var[j])")
                )
      fy11 <- eval(parse(text=pdf11))
pdfArr[,j] <- fy11
if (!is.null(fam$parameters$mu)) 
  {
  m.title <- bquote(paste(.(fname),"(",paste(mu," = ",.(sprintf("%.2f",mu.var[j])),")")))
  }      
#paste0(" mu = ",  sprintf("%.2f",mu.var[j]),    ifelse(!is.null(fam$parameters$sigma),",", " "))
if (!is.null(fam$parameters$sigma)) 
  {
  m.title <- bquote(paste(.(fname),"(",paste(mu," = ",   
                    .(sprintf("%.2f",   mu.var[j])),","),
              paste(sigma," = ",.(sprintf("%.2f",sigma.var[j])),")")))
  }    #      paste(" sigma = ",sprintf("%.2f",sigma.var[j]), ifelse(!is.null(fam$parameters$nu),",", " "), sep="")
if (!is.null(fam$parameters$nu))
  {
  m.title <- bquote(paste(.(fname),"(",paste(mu," = ",   .(sprintf("%.2f",   mu.var[j])),","),
                          paste(sigma," = ",.(sprintf("%.2f",sigma.var[j])),","),
                          paste(nu," = ",   .(sprintf("%.2f",   nu.var[j])),")"),))  
  }
if (!is.null(fam$parameters$tau))
  {
        m.title <- bquote(paste(.(fname),"(",paste(mu," = ",   .(sprintf("%.2f",   mu.var[j])),","),
                                paste(sigma," = ",.(sprintf("%.2f",sigma.var[j])),","),
                                paste(nu, " = ",   .(sprintf("%.2f",   nu.var[j])),","),
                                paste(tau," = ",   .(sprintf("%.2f",   tau.var[j])),")")))    
   }
title.label[j] <- if (!no.title) as.expression(m.title) else ""
if(obj$type=="Discrete")  y.title <- "P(Y=y)" else y.title <- "f(y)"
      plot(y.var , pdfArr[,j], xlab="y",
           ylab=y.title,
           main=title.label[j],
           col=gray(.4),
           frame.plot = frame.plot,
           type=typelh, lty=1, lwd=1, ...)
 }
    par(op)
}###################################################################
else # if only family (NOT A MODEL) is used 
{
      lobs <- max(c(length(mu),length(sigma),length(nu),length(tau)))
if (lobs >= 9)   stop(paste("Use up to eight different combinations of parameters for plotting"))
     plots <- list(c(1,1),c(2,1),c(3,1),c(2,2),c(3,2),c(3,2),c(4,2),c(4,2))
       vvv <- unlist(plots[lobs])
if(any(fname%in%.gamlss.bi.list)) bd <- max
pdfunction <- rep(0, length(y.var))
if(type=="Discrete")   {typelh <- "h" } else   {typelh <- "l" }
if ("mu"%in%names(fam$parameters))
   { if (is.null(mu)) stop("At least one value of mu has to be set")
      mu.var <- rep(mu, length = lobs) 
      if (!fam$mu.valid(mu.var))  stop( "`mu' parameter out of range")
   }
if ("sigma"%in%names(fam$parameters))
   { if (is.null(sigma)) stop("At least one value of sigma has to be set") 
      sigma.var <- rep(sigma, length = lobs)
      if (!fam$sigma.valid(sigma.var))  stop( "`sigma' parameter out of range")
    }
if ("nu"%in%names(fam$parameters))
   { 
      if (is.null(nu)) stop("At least one value of nu has to be set")
      nu.var <- rep(nu, length = lobs)
      if (!fam$nu.valid(nu.var))  stop( "`nu' parameter out of range")
   }
if ("tau"%in%names(fam$parameters))
   { if (is.null(tau)) stop("At least one value of tau has to be set") 
      tau.var <- rep(tau, length = lobs)
      if (!fam$tau.valid(tau.var))  stop( "`tau' parameter out of range")
   }
if (!fam$y.valid(y.var))  stop( "response variable out of range")
     pdfArr <- array(0:0, c(length(y.var),lobs))
title.label <- rep(NA, lobs)     
         op <- par(mfrow=vvv, mar=par("mar")+c(0,0,0,0),
              col.axis=gray(.2),
              col.main=gray(.2),
              col.lab=gray(.2),
              cex=0.6, font=2)
for (j in 1:lobs)
  {
      pdf11 <-  switch (nopar,
                        paste0("d",fname,"(y.var,mu.var[j])"),
                        paste0("d",fname,"(y.var,mu.var[j],sigma.var[j])"),
                        paste0("d",fname,"(y.var,mu.var[j],sigma.var[j], nu.var[j])"),
                        paste0("d",fname,"(y.var,mu.var[j],sigma.var[j], nu.var[j], tau.var[j])")
      )
      fy11 <- eval(parse(text=pdf11))
      pdfArr[,j] <- fy11
      if (!is.null(fam$parameters$mu)) { 
        #bquote(paste(.(fname),"(",paste(mu," = ",.(2),")")))
        m.title <- bquote(paste(.(fname),"(",paste(mu," = ",.(2),")")))} 
      if (!is.null(fam$parameters$sigma)) {
        m.title <-  bquote(paste(.(fname),"(",paste(mu," = ",  .(mu.var[j]),  
                                                    sep=","), 
                                 paste(sigma," = ",.(sigma.var[j])         ),")"))}
      if (!is.null(fam$parameters$nu)) {
        m.title <-  bquote(paste(.(fname),"(",paste(mu," = ", .(mu.var[j]),  
                                                    sep=","), 
                                 paste(sigma," = ",.(sigma.var[j]), sep=","),
                                 paste(   nu," = ",.(   nu.var[j]))        ,")"))}
      if (!is.null(fam$parameters$tau)){
        m.title <-  bquote(paste(.(fname),"(",paste(mu," = ", .(mu.var[j]),  
                                                    sep=","), 
                                 paste(sigma," = ",.(sigma.var[j]), sep=","),
                                 paste(   nu," = ",.(   nu.var[j]), sep=","),
                                 paste(  tau," = ",.(  tau.var[j]))       ,')'))}
      title.label[j] <- if (!no.title) as.expression(m.title) else ""
      
      y.title <- if(type=="Discrete")  "P(Y=y)" else  "f(y)"
      plot(y.var , pdfArr[,j], xlab="y",
           ylab=y.title,
           main= title.label[j],
           col = col,
           frame.plot = frame.plot,
           type = typelh,
           ...)
 }
    par(op)
}
  par(op)
}
