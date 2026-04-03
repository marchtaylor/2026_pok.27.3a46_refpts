
# remove all existing objects ----------------
rm(list = ls())


# required functions ------------------------------------------------------
source("utilities.R")

# data files ----------------------------------------------------------------
load("data/objects.Rdata")

# model settings files ----------------------------------------------------------------
load(file = "model/settings.Rdata")


# redraw stock from SAM uncertainty ---------------------------------------
set.seed(1234) # set seed

# Compute process residuals
samfit0 <- get.procres.sdrep(samfit0, getJointPrecision=TRUE)

# Custom function to generate parameter uncertainty
samfit_unc <- SAM.uncertainty(samfit0, n = it, use_fast_mvrnorm = TRUE)

# expand stock to iters
stock_it <- propagate(stock, it)
dim(stock_it)

for(i in seq(it)){
  stk_i <- FLfse::SAM2FLStock(samfit_unc[[i]], catch_estimate = T)
  iter(stock_it, i) <- stk_i
  print(paste(i, "of", it, "is finished"))
}


# plot(stock_it)


# save --------------------------------------------------------------------

save(stock_it, file = "model/stock_it.Rdata")

