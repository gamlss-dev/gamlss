# ---------------------------------------------------------------------------------------
rqres <- function (pfun = "pNO", 
                   type = c("Continuous", "Discrete", "Mixed"),
               censored = NULL, 
                   ymin = NULL, 
                 mass.p = NULL, 
                prob.mp = NULL,
                      y = y,
                         ... )
{ 
 type <- match.arg(type)
  cdf <- eval(parse(text=pfun))
switch(type, 
  "Continuous"=                   # continous distributions
     {
    rqres <- qnorm(cdf(q=y,...))  # if censored cdf do the right thing
     },  
   "Discrete"=                    # discrete distributions
       {
       # randomize 
        if (is.null(censored))
         {
       aval <- cdf(ifelse(y==ymin,y,y-1), ...) # lower quantile
       aval <- ifelse(y==ymin,0,aval)          # set to zero if y=0
       bval <- cdf(q=y,...)                      # upper quantile
       uval <- runif(length(y),aval,bval) #    gen rand. value
       uval <- ifelse(uval>0.99999,uval-0.1e-10,uval)# 
      rqres <- qnorm(uval)
         }
        else 
         {
        qq <- ifelse(y[,1]==ymin,y[,1], y[,1]-1)    
      aval <- cdf(Surv(qq,y[,2]), ...) # lower quantile
      aval <- ifelse(y[,1]==ymin,0,aval)          # set to zero if y=0
      bval <- cdf(q=y,...)                      # upper quantile
      uval <- runif(length(y[,1]),min=aval,max=bval) #    gen rand. value
      uval <- ifelse(uval>0.99999,uval-0.1e-10,uval)# 
     rqres <- qnorm(ifelse(y[,"status"]==1, uval, bval))
         }
       }, 
   "Mixed"=                     # mixed distributions only two mass points are allowed
       {
      if (is.null(mass.p) && is.null(prob.mp)) 
              stop("For mixed distributions mass.p and prob.mp arguments have to be specified")
      length.mass.p <- length(mass.p)
      #if (length.mass.p != length(prob.mp))
      #   stop("mass.p and prob.mp must have the same length")
      switch(length.mass.p, 
             {                 # it assumes that the fist mass point is at the lower end       
             uval <- ifelse(y==mass.p, runif(length(y),0,prob.mp),cdf(q=y,...))
             },
             {
             uval <- ifelse(y==mass.p[1],runif(length(y),0,prob.mp[,1]),0)
             uval <- ifelse(y>mass.p[1] & y<mass.p[2], cdf(q=y,...),uval)
             uval <- ifelse(y==1,runif(length(y),1-rowSums(prob.mp),1),uval)
             }
             ) 
       rqres <- qnorm(uval)   
       }
       ) 
rqres
}
#----------------------------------------------------------------------------------------
