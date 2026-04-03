# remove all existing objects ----------------
rm(list = ls())


# set seed ----------------------------------------------------------------

set.seed(54321)

# required functions ------------------------------------------------------
source("utilities.R")


# load data and settings ----------------------------------------------------------------

load(file = "data/objects.Rdata", verbose = T)
load(file = "model/settings.Rdata", verbose = T)
load(file = "model/stock_it.Rdata", verbose = T)
load(file = "model/srr_fits_srpars.Rdata", verbose = T)



# BRPs[stk_name,]  
refpts <- FLPar(BRPs["POK",], params = c("Flim", "Fpa", "Fmsy", "lFmsy", "uFmsy", "Fmgt", "Blim", "Bpa", "Btrigger"), units = "NA")




# - CONSTRUCT OM
# GENERATE future deviances: lognormal autocorrelated **
if(use_rhologRec){
  rhologRec <- srpars$rho
}else{
  rhologRec <- srpars$rho
  rhologRec[] <- 0
}

srdevs <- rlnormar1(sdlog = srpars$sigmaR, 
  rho = rhologRec, 
  years = seq(dy, fy), 
  bias.correct = FALSE)

plot(srdevs) + 
  geom_line(data = subset(as.data.frame(srdevs), iter %in% ac(11:15)), 
    mapping = aes(x = year, y = data, group = iter))


# BUILD FLom
tmp <- as.FLSR(stock, model = "segreg")
Bemp <- ssb(tmp)[, rec(tmp) > median(rec(tmp))] |> c() |> sort() |> head(3) |> mean() |> round()

refpts@.Data["Blim",] <- Bemp
om <- FLom(stock = stock_it, refpts = refpts, model = 'mixedsrr',
  params = srpars, deviances = srdevs)

# SETUP om future: average of last 3 years **
om <- fwdWindow(om, end = fy)

# RESAMPLE forecast year slots with Eqsim definitions
# catchability / selectivity
uiter <- seq(it)

if(!sel_const){ # sample years if not const mean value used
  for(i in uiter){
    samp.yrs <- sample(x = sel_years, size = length(iy:fy), replace = T)
    om@stock@harvest[, ac(iy:fy),,,,i] <- om@stock@harvest[,ac(samp.yrs),,,,i]
  }
}else{ # const mean value used
  for(i in uiter){
    mean.yrs <- sel_years
    om@stock@harvest[, ac(iy:fy),,,,i] <- apply(om@stock@harvest[,ac(mean.yrs),,,,i], MARGIN = c(1,3,4,5,6), mean, na.rm = TRUE)
  }
}

# ind weight (stock, landings, discards)
# and biological slots
if(!bio_const){ # sample years if not const mean value used
  for(i in uiter){
    samp.yrs <- sample(x = bio_years, size = length(iy:fy), replace = T)
    
    om@stock@mat[, ac(iy:fy),,,,i] <- om@stock@mat[,ac(samp.yrs),,,,i]
    om@stock@m[, ac(iy:fy),,,,i] <- om@stock@m[,ac(samp.yrs),,,,i]
    om@stock@stock.wt[, ac(iy:fy),,,,i] <- om@stock@stock.wt[,ac(samp.yrs),,,,i]
  
    om@stock@landings.wt[, ac(iy:fy),,,,i] <- om@stock@landings.wt[,ac(samp.yrs),,,,i]
    om@stock@discards.wt[, ac(iy:fy),,,,i] <- om@stock@discards.wt[,ac(samp.yrs),,,,i]
    
    om@stock@landings.n[, ac(iy:fy),,,,i] <- om@stock@landings.n[,ac(samp.yrs),,,,i] / om@stock@catch.n[,ac(samp.yrs),,,,i]
    om@stock@discards.n[, ac(iy:fy),,,,i] <- 1 - om@stock@landings.n[, ac(iy:fy),,,,i]
  }
}else{ # const mean value used
  for(i in uiter){
    mean.yrs <- bio_years

    om@stock@mat[, ac(iy:fy),,,,i] <- apply(om@stock@mat[,ac(mean.yrs),,,,i], MARGIN = c(1,3,4,5,6), mean, na.rm = TRUE)
    om@stock@m[, ac(iy:fy),,,,i] <- apply(om@stock@m[,ac(mean.yrs),,,,i], MARGIN = c(1,3,4,5,6), mean, na.rm = TRUE)

    om@stock@stock.wt[, ac(iy:fy),,,,i] <- apply(om@stock@stock.wt[,ac(mean.yrs),,,,i], MARGIN = c(1,3,4,5,6), mean, na.rm = TRUE)
    om@stock@landings.wt[, ac(iy:fy),,,,i] <- apply(om@stock@landings.wt[,ac(mean.yrs),,,,i], MARGIN = c(1,3,4,5,6), mean, na.rm = TRUE)
    om@stock@discards.wt[, ac(iy:fy),,,,i] <- apply(om@stock@discards.wt[,ac(mean.yrs),,,,i], MARGIN = c(1,3,4,5,6), mean, na.rm = TRUE)
    
    om@stock@landings.n[, ac(iy:fy),,,,i] <- apply(om@stock@landings.n[,ac(mean.yrs),,,,i] / om@stock@catch.n[,ac(mean.yrs),,,,i], MARGIN = c(1,3,4,5,6), mean, na.rm = TRUE)
    om@stock@discards.n[, ac(iy:fy),,,,i] <- 1 - om@stock@landings.n[, ac(iy:fy),,,,i]
  }
}


# SET stochastic rec dy
rec(stock(om))[, ac(dy)] <- rec(om)[1, ac(dy)] * srdevs[, ac(dy)]

# PROJECT forward for iy assumption
ctrl <- fwdControl( 
  data.frame(
    year = c(iy),
    value = c(1),
    quant = c("f"),
    relYear = c(dy)                               
  )
)
om <- fwd(om, control = ctrl)

plot(om)



# F and SSB deviances
# shortcut_devs2 <- function(om, Fcv = 0.212, Fphi = 0.423, SSBcv = 0, SSBphi = 0){
#   devs <- FLQuants(
#     F = rlnormar1(rho = Fphi, years = dimnames(om)$year, n = dims(om)$iter,
#       meanlog = 0, sdlog = Fcv, bias.correct = FALSE),
#     SSB = rlnormar1(rho = 0, years = dimnames(om)$year, n = dims(om)$iter,
#       meanlog = 0, sdlog = SSBcv, bias.correct = FALSE))
#   return(devs)
# }
# sdevs <- shortcut_devs2(om, Fcv = Fcv, Fphi = Fphi, SSBcv = SSBcv)

sdevs <- shortcut_devs(om, Fcv = Fcv, Fphi = Fphi, SSBcv = SSBcv, SSBphi = 0, bias.correct = F)
plot(sdevs$F); mean(sdevs$F)
plot(sdevs$SSB); mean(sdevs$SSB)


# - CONSTRUCT iem, implementation error module not used
iem <- FLiem(method = noise.iem,
  args = list(noise = rlnorm(it, rec(om) %=% 0, 0)))


# - SAVE

save(om, iem, sdevs, file = "model/condition_mse_objects.Rdata")


