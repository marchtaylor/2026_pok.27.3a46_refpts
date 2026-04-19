# remove all existing objects ----------------
rm(list = ls())


# 1. load data ------------------------------------------------------

load(file = "data/objects.Rdata", verbose = T)
load(file = "model/settings.Rdata", verbose = T)
load(file = "model/condition_mse_objects.Rdata", verbose = T)

om@refpts
(Bemp <- om@refpts@.Data["Blim",])
median(om@sr@params@.Data["b",])


# Fscan values
info <- data.frame(Ftarg = Fscan)
info$scen <- seq(nrow(info))
head(info)
dim(info)


# 2. required functions ------------------------------------------------------

# no assessment error, no advice rule (Flim, Btrigger)
parFun_noErr_noAR <- function(x){
  library(mse)
  library(msemodules)
  library(FLasher)
  
  # monitor progress
	sink(file=outfile, append = TRUE)
		cat(paste(x, "of", nrow(infoFmsy), "started @", Sys.time(), "\n"))
	sink()
  
  Ftarg <- infoFmsy$Ftarg[x]
  
  mseargs <- list(iy = iy, fy = fy, data_lag = data_lag, management_lag = management_lag, frq = 1)
  
  # SETUP standard ICES advice rule
  arule <- mpCtrl(list(
    # (est)imation method: shortcut.sa + SSB deviances
    est = mseCtrl(method = shortcut.sa, 
      args = list(SSBdevs = sdevs$SSB*0+1)), # *** SSBcv = 0
    # hcr: hockeystick (fbar ~ ssb | lim, trigger, target, min)
    hcr = mseCtrl(method = hockeystick.hcr, 
      args = list(lim = 0, trigger = 0, target = Ftarg, # *** trigger = 0
        min = 0, metric = "ssb", output = "fbar")),
    # (i)mplementation (sys)tem: tac.is (C ~ F) + F deviances
    isys = mseCtrl(method = tac.is, 
      args = list(recyrs = recyrs, fmin = 0, Fdevs = sdevs$F*0+1)) # *** Fcv = 0 (includes Fphi)
  ))
  
  # - RUN applying ICES advice rule
  system.time(
    advice <- mp(om, iem = iem, ctrl = arule, args = mseargs, parallel = FALSE)
  )
  
  # extract info
  # tracking(advice)
  # inspect(tracking(advice))
  # unique(tracking(advice)$metric)
  res <- list(
    rec = rec(advice@om@stock),
    fbar = fbar(advice@om@stock),
    ssb = ssb(advice@om@stock),
    landings = landings(advice@om@stock),
    catch = catch(advice@om@stock)
  )
  
  # monitor progress
	sink(file = outfile, append = TRUE)
		cat(paste(x, "of", nrow(infoFmsy), "completed @", Sys.time(), "\n"))
	sink()
  
  return(res)
} # end of parFun_noErr_noAR

# with assessment error, no advice rule (Fmsy, Btrigger)
parFun_wErr_noAR <- function(x){
  library(mse)
  library(msemodules)
  library(FLasher)
  
  # monitor progress
	sink(file=outfile, append = TRUE)
		cat(paste(x, "of", nrow(info), "started @", Sys.time(), "\n"))
	sink()
  
  Ftarg <- info$Ftarg[x]
  
  # thus, run MSE, but adapt HCR to Btrigger = 0
  mseargs <- list(iy = iy, fy = fy, data_lag = data_lag, management_lag = management_lag, frq = 1) # management_lag = 0 throws error
  
  # SETUP standard ICES advice rule
  arule <- mpCtrl(list(
    # (est)imation method: shortcut.sa + SSB deviances
    est = mseCtrl(method = shortcut.sa, 
      args = list(SSBdevs = sdevs$SSB*0+1)), # *** SSBcv = 0
    # hcr: hockeystick (fbar ~ ssb | lim, trigger, target, min)
    hcr = mseCtrl(method = hockeystick.hcr, 
      args = list(lim = 0, trigger = 0, target = Ftarg, # *** trigger = 0
        min = 0, metric = "ssb", output = "fbar")),
    # (i)mplementation (sys)tem: tac.is (C ~ F) + F deviances
    isys = mseCtrl(method=tac.is, 
      args = list(recyrs = recyrs, fmin = 0, Fdevs = sdevs$F)) # *** Fcv (includes Fphi)
  ))
  
  # - RUN applying ICES advice rule
  system.time(
    advice <- mp(om, iem = iem, ctrl = arule, args = mseargs, parallel = FALSE)
  )
  
  # extract info
  # tracking(advice)
  res <- list(
    rec = rec(advice@om@stock),
    fbar = fbar(advice@om@stock),
    ssb = ssb(advice@om@stock),
    landings = landings(advice@om@stock),
    catch = catch(advice@om@stock)
  )
  
  # monitor progress
	sink(file = outfile, append = TRUE)
		cat(paste(x, "of", nrow(info), "completed @", Sys.time(), "\n"))
	sink()
  
  return(res)
  
} # end of parFun_wErr_noAR

# with assessment error, with advice rule (Fpa)
parFun_wErr_wAR <- function(x){
  library(mse)
  library(msemodules)
  library(FLasher)
  
  # monitor progress
	sink(file=outfile, append = TRUE)
		cat(paste(x, "of", nrow(info), "started @", Sys.time(), "\n"))
	sink()
  
  Ftarg <- info$Ftarg[x]
  
  # thus, run MSE, but adapt HCR to Btrigger = 0
  mseargs <- list(iy = iy, fy = fy, data_lag = data_lag, management_lag = management_lag, frq = 1) # management_lag = 0 throws error
  
  # SETUP standard ICES advice rule
  arule <- mpCtrl(list(
    # (est)imation method: shortcut.sa + SSB deviances
    est = mseCtrl(method = shortcut.sa, 
      args = list(SSBdevs = sdevs$SSB*0+1)), # *** SSBcv = 0
    # hcr: hockeystick (fbar ~ ssb | lim, trigger, target, min)
    hcr = mseCtrl(method = hockeystick.hcr, 
      args = list(lim = 0, trigger = Btrig_used, target = Ftarg, # *** trigger = Btrigger
        min = 0, metric = "ssb", output = "fbar")),
    # (i)mplementation (sys)tem: tac.is (C ~ F) + F deviances
    isys = mseCtrl(method=tac.is, 
      args = list(recyrs = recyrs, fmin = 0, Fdevs = sdevs$F)) # *** Fcv (includes Fphi)
  ))
  
  # - RUN applying ICES advice rule
  system.time(
    advice <- mp(om, iem = iem, ctrl = arule, args = mseargs, parallel = FALSE)
  )
  
  # extract info
  # tracking(advice)
  res <- list(
    rec = rec(advice@om@stock),
    fbar = fbar(advice@om@stock),
    ssb = ssb(advice@om@stock),
    landings = landings(advice@om@stock),
    catch = catch(advice@om@stock)
  )
  
  # monitor progress
	sink(file = outfile, append = TRUE)
		cat(paste(x, "of", nrow(info), "completed @", Sys.time(), "\n"))
	sink()

  return(res)
} # end of parFun_wErr_wAR


# 3. Run Fscan scenarios --------------------------------------------------

# 3.1. wErr, noAR ------

parFun <- get("parFun_wErr_noAR")

# monitor progress file
outfile = "output_parFun_wErr_noAR.txt"
unlink(outfile)

# only pass relevant ARGs to avoid memory cluttering
ARGS =  c("outfile", "om", "dy", "iy", "fy", "it", "info", "sdevs", "iem", "data_lag", "management_lag", "recyrs")

# Set-up cluster and run
clusterType <- ifelse(Sys.info()["sysname"] == "Windows", "PSOCK", "FORK")
cl <- parallel::makeCluster(ncores, type=clusterType) # set-up cluster
nn <- split(info$scen, info$scen) # indices for each iteration (used as "x" in the parFun)
parallel::clusterExport(cl, varlist = ARGS, envir = environment()) # export info to the cluster
t1 = Sys.time()
RES <- parLapply(cl = cl, X = nn, fun = parFun) # runs the function through each index in nn (index passed to "x" in parFun)
stopCluster(cl) # stops cluster
t2 = Sys.time()
time.taken1 = t2-t1

df1 <- lapply(as.list(seq(info$Ftarg)), function(x){
  tmp1 <- as.data.frame(RES[[x]]$ssb)
  tmp1$slot <- "ssb"
  tmp2 <- as.data.frame(RES[[x]]$landings)
  tmp2$slot <- "landings"
  tmp3 <- as.data.frame(RES[[x]]$catch)
  tmp3$slot <- "catch"
  tmp4 <- as.data.frame(RES[[x]]$rec)
  tmp4$slot <- "rec"
  tmp5 <- as.data.frame(RES[[x]]$fbar)
  tmp5$slot <- "fbar"
  tmp <- rbind(tmp1, tmp2, tmp3, tmp4, tmp5)
  tmp$Ftarg <- info$Ftarg[x]
  return(tmp)
})
df1 <- do.call("rbind", df1)

beepr::beep()



## MSY - wErr, noAR ----
# Fmsy, Bmsy, MSY, FmsyU, FmsyL
aggL <- aggregate(data ~ Ftarg, data = df1, subset = slot == "landings" & year %in% ((fy-ey):fy), FUN = median)
splL <- smooth.spline(x = aggL$Ftarg, y = aggL$data)
newdat <- data.frame(Ftarg = seq(0, max(info$Ftarg), len = 1001))
newdat$landings <- predict(splL, x = newdat$Ftarg)$y
(MSY <- round(newdat$landings[which.max(newdat$landings)]))
(Fmsy <- newdat$Ftarg[which.max(newdat$landings)])
FmsyCI <- range(newdat$Ftarg[which(newdat$landings > MSY*0.95)]) # range of Ftarget +/- Fmsy
hitL <- which(newdat$Ftarg==FmsyCI[1])
hitU <- which(newdat$Ftarg==FmsyCI[2])
(FmsyL <- newdat$Ftarg[hitL])
(FmsyU <- newdat$Ftarg[hitU])

plot(data ~ Ftarg, aggL)
lines(newdat)

# Bmsy
aggB <- aggregate(data ~ Ftarg, data = df1, subset = slot == "ssb" & year %in% ((fy-ey):fy), FUN = median)
spl <- smooth.spline(x = aggB$Ftarg, y = aggB$data)
newdat$B50_wErr_noAR <- predict(spl, x = newdat$Ftarg)$y
(Bmsy <- round(predict(spl, x = Fmsy)$y))
(B0 <- round(subset(aggB, Ftarg == 0)$data))

plot(data ~ Ftarg, aggB, ylab = "SSB")
lines(B50_wErr_noAR ~ Ftarg, newdat)
lines(x = c(-1000, 0, 0), y = c(B0, B0, -1e6), lty = 2, col = 3)
lines(x = c(-1000, Fmsy, Fmsy), y = c(Bmsy, Bmsy, -1e6), lty = 2, col = 4)


# 3.2. noErr, noAR ------

parFun <- get("parFun_noErr_noAR")

infoFmsy <- data.frame(Ftarg = Fmsy)
infoFmsy$scen <- seq(nrow(infoFmsy))

RES <- list()
RES[[1]] <- parFun(x = 1)

df2 <- lapply(as.list(seq(infoFmsy$Ftarg)), function(x){
  tmp1 <- as.data.frame(RES[[x]]$ssb)
  tmp1$slot <- "ssb"
  tmp2 <- as.data.frame(RES[[x]]$landings)
  tmp2$slot <- "landings"
  tmp3 <- as.data.frame(RES[[x]]$catch)
  tmp3$slot <- "catch"
  tmp4 <- as.data.frame(RES[[x]]$rec)
  tmp4$slot <- "rec"
  tmp5 <- as.data.frame(RES[[x]]$fbar)
  tmp5$slot <- "fbar"
  tmp <- rbind(tmp1, tmp2, tmp3, tmp4, tmp5)
  tmp$Ftarg <- info$Ftarg[x]
  return(tmp)
})
df2 <- do.call("rbind", df2)


beepr::beep()


## Btrigger - noErr, noAR ----
aggB05 <- aggregate(data ~ Ftarg, data = df2, subset = slot == "ssb" & year %in% ((fy-ey):fy), FUN = quantile, prob = 0.05)
(Btrig <- round(aggB05$data))

## Bpa - based on Blim. Trumps Btrigger, if larger ----
(Bpa <- round(Bemp * exp(1.645 * 0.2)))

(Btrig_used <- max(Btrig, Bpa)) # used in AR in 3.3


# 3.3. wErr, wAR ------

parFun <- get("parFun_wErr_wAR")

# monitor progress file
outfile = "output_parFun_wErr_wAR.txt"
unlink(outfile)

# only pass relevant ARGs to avoid memory cluttering
ARGS =  c("outfile", "om", "dy", "iy", "fy", "it", "info", "sdevs", "iem", "data_lag", "management_lag", "recyrs", "Btrig_used")

# Set-up cluster and run
clusterType <- ifelse(Sys.info()["sysname"] == "Windows", "PSOCK", "FORK")
cl <- parallel::makeCluster(ncores, type=clusterType) # set-up cluster
nn <- split(info$scen, info$scen) # indices for each iteration (used as "x" in the parFun)
parallel::clusterExport(cl, varlist = ARGS, envir = environment()) # export info to the cluster
t1 = Sys.time()
RES <- parLapply(cl = cl, X = nn, fun = parFun) # runs the function through each index in nn (index passed to "x" in parFun)
stopCluster(cl) # stops cluster
t2 = Sys.time()
time.taken1 = t2-t1

beepr::beep()


df3 <- lapply(as.list(seq(info$Ftarg)), function(x){
  tmp1 <- as.data.frame(RES[[x]]$ssb)
  tmp1$slot <- "ssb"
  tmp2 <- as.data.frame(RES[[x]]$landings)
  tmp2$slot <- "landings"
  tmp3 <- as.data.frame(RES[[x]]$catch)
  tmp3$slot <- "catch"
  tmp4 <- as.data.frame(RES[[x]]$rec)
  tmp4$slot <- "rec"
  tmp5 <- as.data.frame(RES[[x]]$fbar)
  tmp5$slot <- "fbar"
  tmp <- rbind(tmp1, tmp2, tmp3, tmp4, tmp5)
  tmp$Ftarg <- info$Ftarg[x]
  return(tmp)
})
df3 <- do.call("rbind", df3)

df3 <- tibble(df3)



## Fpa (uses iter-specific empirical Blim) ----
df_blim <- data.frame(blim = c(om@sr@params["b"]), iter = ac(seq(c(om@sr@params["b"]))))
df3a <- df3 |> left_join(df_blim)
df3a <- df3a |> mutate(unsafe = if_else(slot == "ssb", data < blim, NA))

tmp <- df3a |> filter(slot == "ssb", year %in% ((fy-ey):fy)) |> group_by(Ftarg, year) |> 
  summarise(annP = mean(data < blim)) |> ungroup()
tmp <- tmp |> group_by(Ftarg) |> summarise(p1 = mean(annP), p3 = max(annP))
tmp <- tmp |> pivot_longer(cols = -1, names_to = "type", values_to = "risk")
tmp

newdat <- data.frame(Ftarg = seq(0,max(df3$Ftarg), by = 0.001), type = "p3")
newdat$risk <- with(subset(tmp, type == "p3"), approx(x = Ftarg, y = risk, xout = newdat$Ftarg, rule = 2))$y
(Fpa <- newdat |> group_by(type) |> summarise(Fpa = max(Ftarg[risk < 0.05])) |> pull(Fpa))

ggplot(tmp) + aes(x = Ftarg, y = risk, group = type, colour = type) + 
  geom_hline(yintercept = 0.05, lty = 3) + 
  geom_vline(xintercept = Fpa, lty = 3) + 
  geom_line() +
  geom_point() +
  geom_line(data = newdat, lty = 2, color = 4)





# aggB05 <- aggregate(data ~ Ftarg, data = df3, subset = slot == "ssb" & year %in% ((fy-ey):fy), FUN = quantile, prob = 0.05)
# # interpol.F <- approxfun(x = aggB05$Ftarg, y = aggB05$data)
# # newdat$B05_wAR <- interpol.F(v = newdat$Ftarg)
# spl <- smooth.spline(x = aggB05$Ftarg, y = aggB05$data)
# newdat$B05_noErr_wAR <- predict(spl, x = newdat$Ftarg)$y
# (Fpa <- max(newdat$Ftarg[newdat$B05_noErr_wAR > Bemp]))
# 
# plot(B05_noErr_wAR ~ Ftarg, newdat, t = "l")
# abline(h = Bemp, lty = 2)
# abline(v = Fpa, lty = 2)




# summary BRPs ------------------------------------------------------------
newBRPs <- BRPs*NA
# newBRPs[stk_name, "Flim"] <- Flim
newBRPs[, "Fpa"] <- Fpa
newBRPs[, "Fmsy"] <- Fmsy
newBRPs[, "FmsyL"] <- FmsyL
newBRPs[, "FmsyH"] <- FmsyU
newBRPs[, "Btrigger"] <- Btrig
newBRPs[, "Blim"] <- round(Bemp)
newBRPs[, "Bpa"] <- Bpa
newBRPs[, "MSY"] <- MSY
newBRPs[, "Bmsy"] <- Bmsy
newBRPs[, "B0"] <- B0
newBRPs



# save results ----------------------------------------------------------


save(newBRPs, file = "model/newBRPs.Rdata")
saveRDS(df1, file = "model/sum_wErr_noAR.rds")
saveRDS(df2, file = "model/sum_noErr_noAR.rds")
saveRDS(df3, file = "model/sum_wErr_wAR.rds")

# load(file = "model/newBRPs.Rdata")
# df1 <- readRDS(file = "model/sum_wErr_noAR.rds")
# df2 <- readRDS(file = "model/sum_noErr_noAR.rds")
# df3 <- readRDS(file = "model/sum_wErr_wAR.rds")

