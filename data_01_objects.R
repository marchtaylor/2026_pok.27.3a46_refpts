

# remove all existing objects ----------------
rm(list = ls())



load(file = "boot/data/BRPs.Rdata") # old reference points
load(file = "boot/data/eqsim_lut.Rdata")  # old eqsim settings

# load stock object
tmp <- load("boot/data/2025_update_fit.Rdata", verbose = T)
samfit0 <- get(tmp); rm(fit)

stock <- FLfse::SAM2FLStock(samfit0, catch_estimate = T)
units(stock) <- standardUnits(stock)
DIMS <- dims(stock)
RANGE <- range(stock)



# save  -------------------------------------------------------------------

save(stock, samfit0, eqsim_lut, BRPs, file = "data/objects.Rdata")

