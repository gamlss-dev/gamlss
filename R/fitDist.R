# MS and VV 25-9-11
# Fitting a group of gamlss.family distributions
# TO DO
# add all available distributions
#-------------------------------------------------------------------------------
# the grouping of distributions
# #-------------------------------------------------------------------------------
# .realline <- c( "NO", "GU", "RG" ,"LO", "NET",     # 2 par
#   "TF", "TF2", "PE","PE2", "SN1", "SN2", "exGAUS", # 3 par
#   "SHASH", "SHASHo","SHASHo2",                     # 4 par
#   "EGB2", "JSU", "JSUo",                           # 4 par 
#   "SEP1", "SEP2", "SEP3", "SEP4",                  # 4 par 
#    "ST1", "ST2", "ST3", "ST4", "ST5", "SST",       # 4 par 
#   "GT")  
# #--------------------------------------------------------------------------------
# .realplus <- c("EXP",                            # 1 par
#   "GA","IG","LOGNO", "WEI3", "IGAMMA","PARETO2", # 2 par
#   "BCCGo", "exGAUS", "GG", "GIG",                # 3 par
#   "BCTo", "BCPEo", "GB2")                        # 4 par 
# #--------------------------------------------------------------------------------
# .real0to1 <- c("BE", "LOGITNO", "BEINF", "BEINF0", "BEINF1", "BEOI", "BEZI", "GB1")
# #--------------------------------------------------------------------------------
# .realAll <- union(.realline, .realplus)
# # .realAll <- c( "EXP",                                      # 1 par
# #   "GA","IG","LOGNO", "WEI3", "IGAMMA","PARETO2","GP",      # 2 par
# #   "BCCGo", "exGAUS", "GG", "GIG",                          # 3 par
# #   "BCTo", "BCPEo",                                         # 4 par
# #   "GU", "RG" ,"LO", "NET", # 2 par
# #   "TF", "PE", "SN1", "SN2",    # 3 par
# #   "SHASHo", "EGB2", "JSU", "SEP1", "SEP2", "SEP3", "SEP4", # 4 par
# #    "ST1", "ST2", "ST3", "ST4", "ST5", "GT" )               # 4 par
# #           
# 
# #--------------------------------------------------------------------------------               
# #--------------------------------------------------------------------------------
# .counts <- c("GEOM", "GEOMo", "LG", "PO", "YULE","ZIPF", # 1 par
#              "NBI", "NBII", "PIG",  "WARING", "ZALG", 
#              "ZAP", "ZAZIPF", "ZIP", "ZIP2", "GPO", "DPO",# 2 par 
#              "BNB", "NBF", "DEL", "SICHEL", "SI",
#              "ZANBI", "ZAPIG", "ZINBI", "ZIPIG", "ZANBF", # 3 par
#              "ZABNB", "ZASICHEL", "ZINBF", "ZIBNB",       # 4 par
#              "ZISICHEL")
# .counts <- c("PO", "GEOM", "GEOMo","LG", "YULE", "ZIPF", 
#               "WARING", "GPO", "DPO", "BNB", "NBF",      # 
#               "NBI", "NBII", "PIG", "ZIP","ZIP2", "ZAP", "ZALG", # 2 par
#               "DEL", "ZAZIPF",
#               "SI", "SICHEL","ZANBI",  "ZAPIG", "ZINBI",  "ZIPIG", "ZINBF",
#               "ZABNB", "ZASICHEL"
# #--------------------------------------------------------------------------------
# .binom <- c("BI",                            # 1 par
#             "BB", "ZIBI", "ZABI", "DBI",     # 2 par
#             "ZIBB",  "ZABB" )                # 3 par
#-------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------
fitDist <- function(y,
                    k = 2, # for the AIC
                 type = c("realAll", "realline", "realplus","real0to1","counts", "binom" ), 
           try.gamlss = FALSE,  # whether to try the gamlss() if gamlssML() fails
                extra = NULL,  # for extra distributions to include 
                 data = NULL, trace = FALSE, ...)
{
  # if (!is.null(data)) {attach(data); on.exit(detach(data))}
  #if (!is.null(data)) {attach(data, name="TheDatA"); on.exit(detach(TheDatA))}
       y <- if (!is.null(data)) get(deparse(substitute(y)), envir=as.environment(data)) else y
    type <- match.arg(type)
    DIST <- switch(type, "realAll"=.realAll, 
                        "realline"=.realline, 
                        "realplus"=.realplus,
                        "real0to1"=.real0to1,
                          "counts"=.counts,
                           "binom"=.binom 
                  )
if  (!is.null(extra)) DIST <- unique(c(DIST, extra))
    # we need weights here 
#if  ("weights"%in%names(list(...))) wlist(...)$weights else rep()
    m0 <- switch(type,  "realAll"= gamlssML(y, family=NO, ...),
                       "realline"= gamlssML(y, family=NO, ...), 
                       "realplus"= gamlssML(y, family=EXP, ...),
                       "real0to1"= gamlssML(y, family=BE, ...),
                         "counts"= gamlssML(y, family=PO, ...),
                          "binom"= gamlssML(y, family=BI, ...) 
                 ) 
  failed <- list() 
    fits <- list()
#      ow <- options("warn")
if (trace)     cat("---------------------------------------- ","\n")
if (trace)     cat("fitting different", type, "distributions", "\n")
    pb <- txtProgressBar(max = length(DIST), style=3)
    for (i in 1:length(DIST)) 
{
      setTxtProgressBar(pb, i)  
    m1 <- try(gamlssML(y,family=DIST[i], ...), silent=TRUE)
        if (any(class(m1)%in%"try-error")&&try.gamlss==TRUE) 
        { 
         m1 <-  try(gamlss(y~1,family=DIST[i], trace=FALSE, ...),  silent=TRUE)
        }
    
       if (any(class(m1)%in%"try-error"))
       {
            failed <- c(failed, DIST[i]) 
            if (trace)     cat(i, " ", DIST[i], "FAILED", "\n")
       }
       else
      {
               aic <- AIC(m1, k=k)
        names(aic) <- DIST[i]
        if (trace)     cat(i, " ", DIST[i], aic, "\n")
       # options(warn = 1)  
        fits <- c(fits, aic)
        if (AIC(m1, k=k) < AIC(m0, k=k)) 
        {
          m0<-m1 
        }
      }
}
    close(pb)    
#options(ow)     
 m0$failed <- failed
      fits <- unlist(fits)
   m0$fits <- fits[order(fits)]
m0  
}
#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
