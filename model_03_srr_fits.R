# remove all existing objects ----------------
rm(list = ls())

# data files ----------------------------------------------------------------
load("data/objects.Rdata", verbose = T)

# model settings files ----------------------------------------------------------------
load(file = "model/settings.Rdata", verbose = T)

# propagated stock file ----------------------------------------------------------------
load(file = "model/stock_it.Rdata", verbose = T)





# stock-recruitment relationship(s) ------------------------------------

SRmods <- srr_models #c("bevholt", "ricker", "segreg3")
SRmods <- tolower(SRmods)
SRmodsAdj <- replace(SRmods, SRmods == "segreg3", "segreg")

# truncate to srr years
stkSRyrs <- stock[,ac(srr_years)]

# Define single B_empirical (Bemp)
tmp <- as.FLSR(stock, model = "segreg")
Bemp <-  ssb(tmp)[, rec(tmp) > median(rec(tmp))] |> c() |> sort() |> head(3) |> mean()

# fit one iteration
if("segreg3" %in% SRmods){
  srpars <- bootstrapSR(stkSRyrs, iters = 1,
    models = SRmodsAdj, method = "best", Blim = Bemp,
    s = 0.7, s.logitsd = 0.4)
}else{
  srpars <- bootstrapSR(stkSRyrs, iters = 1,
    models = SRmodsAdj, method = "best",
    s = 0.7, s.logitsd = 0.4)
}


# expand srpars to iters
dim(srpars)
srpars <- propagate(srpars, it)
dim(srpars)


# fit to each stock iteration (FLSRTMB) ---------------------------------------------------

# for(i in seq(it)){
#   # stk_i <- as.FLStock.sam(samfit_unc[[i]])
#   stk_i <- iter(stock_it, i)
#   sr_i <- as.FLSR(stk_i, model = SRmodsAdj) # temporary for calc of Bemp_i
# 
#   # Type 5 Blim: average of lowest 3 SSB resulting in R > median(R)
#   # KEY: Taken from entire time series (WKBGAD)
#   Bemp_i <-  ssb(sr_i)[, rec(sr_i) > median(rec(sr_i))] |> c() |> sort() |> head(3) |> mean()
# 
#   # Truncate historical time series for fitting SR
#   trunc_stk <- stk_i[,ac(srr_years)]
#   
#   # Fit SR
#   if("segreg3" %in% SRmods){
#     srpars_i <- bootstrapSR(trunc_stk, iters = 1,
#       models = SRmodsAdj, method = "best", Blim = Bemp_i,
#       s = 0.7, s.logitsd = 0.4)
#   }else{
#     srpars_i <- bootstrapSR(trunc_stk, iters = 1,
#       models = SRmodsAdj, method = "best",
#       s = 0.7, s.logitsd = 0.4)
#   }
#   
#   # pass pars to srpars
#   srpars@.Data[,i] <- srpars_i
#   
#   print(paste(i, "of", it, "is finished"))
# }


parFun <- function(i){
  library(FLCore)
  library(FLSRTMB)
  
  # stk_i <- as.FLStock.sam(samfit_unc[[i]])
  stk_i <- iter(stock_it, i)
  sr_i <- as.FLSR(stk_i, model = SRmodsAdj) # temporary for calc of Bemp_i

  # Type 5 Blim: average of lowest 3 SSB resulting in R > median(R)
  # KEY: Taken from entire time series (WKBGAD)
  Bemp_i <-  ssb(sr_i)[, rec(sr_i) > median(rec(sr_i))] |> c() |> sort() |> head(3) |> mean()

  # Truncate historical time series for fitting SR
  trunc_stk <- stk_i[,ac(srr_years)]
  
  # Fit SR
  if("segreg3" %in% SRmods){
    srpars_i <- bootstrapSR(trunc_stk, iters = 1,
      models = SRmodsAdj, method = "best", Blim = Bemp_i,
      s = 0.7, s.logitsd = 0.4)
  }else{
    srpars_i <- bootstrapSR(trunc_stk, iters = 1,
      models = SRmodsAdj, method = "best",
      s = 0.7, s.logitsd = 0.4)
  }
  
  # return
  return(srpars_i)
}


# only pass relevant ARGs to avoid memory cluttering
ARGS =  c("stock_it", "SRmods", "SRmodsAdj", "srr_years")

# Set-up cluster and run
clusterType <- ifelse(Sys.info()["sysname"] == "Windows", "PSOCK", "FORK")
cl <- parallel::makeCluster(ncores, type=clusterType) # set-up cluster
nn <- split(seq(it), seq(it)) # indices for each iteration (used as "x" in the parFun)
parallel::clusterExport(cl, varlist = ARGS, envir = environment()) # export info to the cluster
t1 = Sys.time()
RES <- parLapply(cl = cl, X = nn, fun = parFun) # runs the function through each index in nn (index passed to "x" in parFun)
stopCluster(cl) # stops cluster
t2 = Sys.time()
time.taken1 = t2-t1

# put results in srpars
for(i in seq(it)){
  srpars@.Data[,i] <- RES[[i]]
}

class(srpars)
# plot(srpars)
srpars_redraw_FLSRTMB <- srpars



# save --------------------------------------------------------------------

save(srpars, file = "model/srr_fits_srpars.Rdata")
  



