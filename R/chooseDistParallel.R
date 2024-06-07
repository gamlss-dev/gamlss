# new function to fit different distributions
# Authors Mikis Stasinopoulos and Vlasios Voudouris
# TO DO
# i)  should allow for multiple k (OK)
# ii) need paralell programing (Partly Ok it works 
#      but I am not sure what happends when use  widows)
# iii) cat() can be better (OK) but it does not works for parallel
# vi) output should be a matrix (OK)) with some functionality (see odrered function) 
# v)  create new list for fitting all possible distribution with different 
#    parametrizations (OK)
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
################################################################################
# this grouping was checked on the 27-4-18
#-------------------------------------------------------------------------------
# the grouping of distributions
#------------------------------------------------------------------------------- 
# this
# group of distribution with interval ranging from -infinity to +infinity
.realline <- c( "NO", "GU", "RG" ,"LO", "NET",     # 2 par
                    "TF", "TF2", "PE","PE2", "SN1", "SN2", "exGAUS", # 3 par
                    "SHASH", "SHASHo","SHASHo2",                     # 4 par
                    "EGB2", "JSU", "JSUo",                           # 4 par 
                    "SEP1", "SEP2", "SEP3", "SEP4",                  # 4 par 
                    "ST1", "ST2", "ST3", "ST4", "ST5", "SST",        # 4 par 
                     "GT")  
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# Group of distribution with interval ranging from 0  to +infinity
.realplus <- c( "EXP", # 1 par
                "GA","IG","LOGNO", "LOGNO2","WEI", "WEI2", "WEI3", "IGAMMA",
                "PARETO2", "PARETO2o", "GP", # 2 par
                "BCCG", "BCCGo", "exGAUS", "GG", "GIG", "LNO",  # 3 par
                "BCTo", "BCT", "BCPEo", "BCPE", "GB2")  # 4 par 
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# Group of distribution with interval ranging from 0 to 1
.real0to1 <- c("BE", "BEo", # 2 par
               "BEINF0", "BEINF1", "LOGITNO", "SIMPLEX", #2 par 
               "BEOI", "BEZI", # 3 par
               "BEINF", # 4 par
               "GB1") # par

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# Group of distribution with interval ranging from -infinity to +infinity (.realline) or 0 to +infinity (.realplus)
.realAll <- union(.realline, .realplus)
# .realAllALL <- c( "EXP", # 1 par
#                "GA","IG","LOGNO", "LOGNO2", "WEI", "WEI2", "WEI3", "IGAMMA",
#                "PARETO2", "PARETO2o", "GP", # 2 par
#                "BCCG", "BCCGo",  "GG", "GIG", "LNO",  # 3 par
#                "BCTo", "BCT", "BCPEo", "BCPE", "GB2",
#                "NO","GU", "RG" ,"LO", "NET", # 2 par
#                "TF", "TF2", "PE", "PE2", "SN1", "SN2", "exGAUS",   # 3 par
#                "SHASH", "SHASHo","SHASHo2", "EGB2", "JSU", "JSUo", 
#                "SEP", "SEP1", "SEP2", "SEP3", "SEP4", "SEP", # 4 par
#                "ST1", "ST2", "ST3", "ST3C", "ST4", "ST5", "SST", "GT")  # 4 par

#-------------------------------------------------------------------------------               
#-------------------------------------------------------------------------------
# Group of distribution for counting
.counts <- c("PO", "GEOM", "GEOMo","LG", "YULE", "ZIPF", # 1 par
             "WARING", "GPO", "DPO", "BNB", "NBF",       # 
             "NBI", "NBII", "PIG", "ZIP","ZIP2", "ZAP", "ZALG", # 2 par
             "DEL", "ZAZIPF", "SI", "SICHEL","ZANBI",  "ZAPIG", 
             "ZINBI",  "ZIPIG", "ZINBF",
             "ZABNB", "ZASICHEL", "ZINBF",  "ZIBNB", "ZISICHEL")
             
#--------------------------------------------------------------------------------
# Binomial group of distributions
.binom <- c("BI",  # 1 par
            "BB", "DBI", "ZIBI", "ZABI", # 2 par
            "ZIBB", "ZABB")
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
chooseDist <- function(object,
                       k = c(2, 3.84, round(log(length(object$y)),2)), # for the AIC
                    type = c("realAll", "realline", "realplus","real0to1","counts", "binom", "extra" ), 
                   extra = NULL,  # for extra distributions to include 
                   trace = FALSE,
                parallel = c("no", "multicore", "snow"), #
                   ncpus = 1L, #integer: number of processes to be used in parallel operation: typically one would chose this to the number of available CPUs
                      cl = NULL, # An optional parallel or snow cluster for use if parallel = "snow". If not supplied, a cluster on the local machine is created for the duration of the boot call.
                       ...)
{
  ## get the type of distribution  
  type <- match.arg(type)
  DIST <- switch(type, "realAll"= .realAll, 
                      "realline"= .realline, 
                      "realplus"= .realplus,
                      "real0to1"= .real0to1,
                        "counts"= .counts,
                         "binom"= .binom, 
                         "extra"= extra)
  if (type=="extra"&&is.null(extra)) stop("extra is not set")
  if  (!is.null(extra)) DIST <- unique(c(DIST, extra))
  ##   
       m0 <- object
  klength <- length(k)
  AiC  <- rep(NA, klength) #matrix(NA, nrow=length(DIST), ncol= klength, dimnames=list(DIST,  as.character(k)))
#--------------- PARALLEL-------------------------------------------------------
#----------------SET UP PART----------------------------------------------------
    parallel <- match.arg(parallel)
     have_mc <- have_snow <- FALSE
    if (parallel != "no" && ncpus > 1L) 
      {
       if (parallel == "multicore") 
         have_mc <- .Platform$OS.type != "windows"
       else if (parallel == "snow") 
         have_snow <- TRUE
       if (!have_mc && !have_snow) 
         ncpus <- 1L
       loadNamespace("parallel")
     } 
# -------------- finish parallel------------------------------------------------
# define the function 
fun <- function(dist,...)
  {    
    m1 <- try(update(object,family=dist, trace=FALSE,...), silent=TRUE)
    if (any(class(m1)%in%"try-error")) 
    { 
    m1 <- try(update(object,family=dist, trace=FALSE, ...),  silent=TRUE)
    }
    else
    {
      for (j in 1:klength) AiC[j] <- AIC(m1, k=k[j])
      if (trace)     cat(dist, "\n",AiC, "\n")
    }
    c(AiC)        # autput of the function
  }
#----------------------------------------------------------------    
#----------------------------------------------------------------
#----------------------------------------------------------------
# --------  parallel --------------------------------------------
  MM <- if (ncpus > 1L && (have_mc || have_snow)) 
     {
       if (have_mc) 
         {# sapply(scope, fn)
         matrix(unlist(parallel::mclapply(DIST, fun, mc.cores = ncpus)), 
                ncol=klength, byrow = T, dimnames = list(DIST,  as.character(k)) )   
       }
      else if (have_snow) 
        {
           list(...)
          if (is.null(cl)) 
            { # make the cluster
            if (.Platform$OS.type == "windows")
            {
              cl <- parallel::makePSOCKcluster(rep("localhost", ncpus))
              clusterEvalQ(cl,pacman::p_load(gamlss)) 
              exp.data =  paste0(object$call$data)
              clusterExport(cl, c(ls(envir = .GlobalEnv), exp.data))
            } else  cl <- parallel::makeForkCluster(ncpus)
             if (RNGkind()[1L] == "L'Ecuyer-CMRG") 
                 parallel::clusterSetRNGStream(cl)
             res <-  matrix(unlist((parallel::parLapply(cl, DIST, fun))),
                                  ncol=klength, byrow = T, 
                            dimnames = list(DIST,  as.character(k))) 
             parallel::stopCluster(cl)
             res
            } 
        else parallel::parLapply(cl, DIST, fun)# use existing cluster
        }
      }# end parallel -----
     else  matrix(sapply(DIST, fun, ...), ncol=klength, byrow = T, 
                    dimnames = list(DIST,  as.character(k)))  
#----------------------------------------------------------------  
#----------------------------------------------------------------    
     for (i in 1:length(k)) cat("minimum GAIC(k=",k[i],") family:", 
                                names(which.min(MM[,i])), "\n")
#----------------------------------------------------------------
     MM  
}
#----------------------------------------------------------------
getOrder <- function(obj, column=1) 
{
  if (!is.matrix(obj)) stop("the object should be a matrix")
  out <- obj[,column][order(obj[,column])]
  name <- colnames(obj)[column]
  cat("GAIG with k=", name, "\n" )
  out 
}
 
#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
# chooseDistPred <- function(object,
#                            newdata = NULL,  
#                               rand = NULL,
#                               type = c("realAll", "realline", "realplus","real0to1","counts", "binom", "extra" ), 
#                              extra = NULL,  # for extra distributions to include 
#                              trace = FALSE, ...)
# {
#   ## get the type of distribution 
#       newData <- if(is.null(newdata)) FALSE else TRUE 
#          type <- match.arg(type)
#   if (is.null(rand)&&is.null(newdata)) stop("rand or newdata should be set")
#   if (!is.null(rand)&&!is.null(newdata)) stop("only rand or newdata should be set NOT both")
#   DIST <- switch(type, "realAll"=.realAll, 
#                       "realline"=.realline, 
#                       "realplus"=.realplus,
#                       "real0to1"=.real0to1,
#                         "counts"=.counts,
#                          "binom"=.binom, 
#                          "extra"=extra)
#   if (type=="extra"&&is.null(extra)) stop("extra is not set")
#   if  (!is.null(extra)) DIST <- unique(c(DIST, extra))
#   ##   
#        m0 <- object
#      tgd0 <- getTGD(m0, newdata=newdata)
#      AiC  <- NA #matrix(NA, nrow=length(DIST), ncol= klength, dimnames=list(DIST,  as.character(k)))
#   # define the function 
# 
#   fun <- function(dist)
#   {    
#     m1 <- try(update(object,family=dist, trace=FALSE,...), silent=TRUE)
#     if (any(class(m1)%in%"try-error")) 
#      { 
#       aic <-   NA
#       }
#     else
#     {
#       tgd1 <- getTGD(m1, newdata=newdata)
#       if (trace)     cat(dist, "\n", tgd1$TGD, "\n")
#       aic <- tgd1$TGD    
#       if (TGD(tgd1) < TGD(tgd0)) 
#       {
#         m0 <<- m1 # saving the best model according to k[order.by]
#       }
#     }
#     aic        # autput of the function
#   }
#   #----------------------------------------------------------------    
#   #  should be parallise
#   MM <- sapply(DIST, fun) 
#   #----------------------------------------------------------------  
#   ## save it in the final model   
#   m0$TGD <- MM[order(MM)]
#   #----------------------------------------------------------------    
#   #----------------------------------------------------------------    
#   return(m0)  
# }
# #--------------------------------------------------------------------------------------
# #--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
chooseDistPred <- function(object,
                           type = c("realAll", "realline", "realplus","real0to1","counts", "binom", "extra" ), 
                           extra = NULL,  # for extra distributions to include 
                           trace = FALSE,
                           parallel = c("no", "multicore", "snow"), #
                           ncpus = 1L, #integer: number of processes to be used in parallel operation: typically one would chose this to the number of available CPUs
                           cl = NULL, # An optional parallel or snow cluster for use if parallel = "snow". If not supplied, a cluster on t     
                           newdata = NULL,
                           rand = NULL,
                           ...)
{
  ## get the type of distribution
  newData <- if(is.null(newdata)) FALSE else TRUE
  type <- match.arg(type)
  if (is.null(rand)&&is.null(newdata)) stop("rand or newdata should be set")
  if (!is.null(rand)&&!is.null(newdata)) stop("only rand or newdata should be set NOT both")
  DIST <- switch(type, "realAll"=.realAll,
                 "realline"=.realline,
                 "realplus"=.realplus,
                 "real0to1"=.real0to1,
                 "counts"=.counts,
                 "binom"=.binom,
                 "extra"= extra)
  if (type=="extra"&&is.null(extra)) stop("extra is not set")
  if  (!is.null(extra)) DIST <- unique(c(DIST, extra))
  ##
  m0 <- object
  tgd0 <- getTGD(m0, newdata=newdata)
  AiC  <- rep(NA, 1)
  #--------------- PARALLEL-------------------------------------------------------
  #----------------SET UP PART----------------------------------------------------
  parallel <- match.arg(parallel)
  have_mc <- have_snow <- FALSE
  if (parallel != "no" && ncpus > 1L) 
  {
    if (parallel == "multicore") 
      have_mc <- .Platform$OS.type != "windows"
    else if (parallel == "snow") 
      have_snow <- TRUE
    if (!have_mc && !have_snow) 
      ncpus <- 1L
    loadNamespace("parallel")
  } 
  # -------------- finish parallel------------------------------------------------     
  # define the function
  fun <- function(dist)
  {
    m1 <- try(update(object,family=dist, trace=FALSE,...), silent=TRUE)
    if (any(class(m1)%in%"try-error"))
    {
      m1 <- try(update(object,family=dist, trace=FALSE, ...),  silent=TRUE)
    }
    else
    {
      tgd1 <- getTGD(m1, newdata=newdata)
      if (trace)     cat(dist, "\n", tgd1$TGD, "\n")
      AiC <- tgd1$TGD
      # if (TGD(tgd1) < TGD(tgd0))
      # {
      #   m0 <<- m1 # saving the best model according to k[order.by]
      # }
    }
    c(AiC)        # autput of the function
  }
  #----------------------------------------------------------------    
  #----------------------------------------------------------------
  #----------------------------------------------------------------
  # --------  parallel --------------------------------------------
  MM <- if (ncpus > 1L && (have_mc || have_snow)) 
  {
    if (have_mc) 
    {# sapply(scope, fn)
      unlist(parallel::mclapply(DIST, fun, mc.cores = ncpus))
    }
    else if (have_snow) 
    {
      list(...)
      if (is.null(cl)) 
      { # make the cluster
        if (.Platform$OS.type == "windows")
        {
          cl <- parallel::makePSOCKcluster(rep("localhost", ncpus))
          clusterEvalQ(cl,pacman::p_load(gamlss)) 
          exp.data =  paste0(object$call$data)
          clusterExport(cl, c(ls(envir = .GlobalEnv), exp.data))
        } else cl <- parallel::makeForkCluster(ncpus)
        if (RNGkind()[1L] == "L'Ecuyer-CMRG") 
          parallel::clusterSetRNGStream(cl)
        res <-  unlist((parallel::parLapply(cl, DIST, fun)))
        parallel::stopCluster(cl)
        res
      } 
      else parallel::parLapply(cl, DIST, fun)# use existing cluster
    }
  }# end parallel -----
  else  sapply(DIST, fun) 
  #----------------------------------------------------------------  
  #----------------------------------------------------------------         
  names(MM) <- DIST
  #----------------------------------------------------------------
  ## save it in the final model
  #m0$TGD <- MM[order(MM)]
  #----------------------------------------------------------------
  #----------------------------------------------------------------
  MM[order(MM)]
}
#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------

