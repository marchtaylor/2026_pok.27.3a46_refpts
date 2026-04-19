# remove all existing objects ----------------
rm(list = ls())



# 1. load data ------------------------------------------------------

load(file = "data/objects.Rdata", verbose = T)
load(file = "model/settings.Rdata", verbose = T)
load(file = "model/condition_mse_objects.Rdata", verbose = T)
load(file = "model/srr_fits_srpars.Rdata", verbose = T)


# 2. required functions ------------------------------------------------------

segreg3 <- function(ab, ssb) log(ifelse(ssb >= Bemp, ab$a * Bemp, ab$a * ssb))


# 3. other settings -------------------------------------------------------------

# Define single B_empirical (Bemp)
tmp <- as.FLSR(stock, model = "segreg")
Bemp <-  ssb(tmp)[, rec(tmp) > median(rec(tmp))] |> c() |> sort() |> head(3) |> mean() |> round()
Bpa <- round(Bemp * exp(1.645 * 0.2))

remove.years <- an(dimnames(stock)$year)[!an(dimnames(stock)$year) %in% srr_years]
 


# 4. Eqsim with default ------------------------------------------------

# fit srrs (bootstrapping with msy package)
srr_fit <- eqsr_fit(stk = stock, nsamp = it, 
  models = c("segreg3"), remove.years = remove.years) 


# # truncate to srr years
# stkSRyrs <- stock[,ac(srr_years)]
# 
# srpar_tmp <- bootstrapSR(stkSRyrs, iters = 1,
#     models = "segreg", Blim = Bemp)
# 
# c(srpar_tmp["rho"])
# 
# use_FLSR_pars <- FALSE
# if(use_FLSR_pars){
# }


## 4.1. wErr_noAR -------------------------------


set.seed(1111)
FIT1 <- eqsim_run(srr_fit, 
  bio.years = range(bio_years), 
  sel.years = range(sel_years),
  bio.const = bio_const, 
  sel.const = sel_const,
  rhologRec = mean(srpars$rho),
  Fcv = Fcv, Fphi = Fphi, SSBcv = SSBcv,
  Btrigger = 0, Blim = Bemp, Bpa = Bpa,
  Fscan = Fscan,
  verbose = FALSE,
  # Nrun = c(fy-iy), 
  keep.sims = TRUE
)



### MSY ref points ----
# Fmsy, Bmsy, MSY, FmsyU, FmsyL

(Fmsy <- round(FIT1$Refs2["lanF","medianMSY"], 3))
(FmsyL <- round(FIT1$Refs2["lanF","Medlower"], 3))
(FmsyU <- round(FIT1$Refs2["lanF","Medupper"], 3))
(Bmsy <- round(FIT1$Refs2["lanB","medianMSY"]))
(MSY <- round(FIT1$Refs2["landings","medianMSY"]))
(B0 <- (FIT1$rbp |> filter(variable == "Spawning stock biomass", Ftarget == 0) |> select(p50))[,1])

eqsim_plot_range(FIT1, type="median")



## 4.2. noErr_noAR -----------------------------

set.seed(1111)
FIT2 <- eqsim_run(srr_fit, 
  bio.years = range(bio_years), 
  sel.years = range(sel_years),
  bio.const = bio_const, 
  sel.const = sel_const,
  rhologRec = mean(srpars$rho),
  Fcv = 0, Fphi = 0, SSBcv = 0,
  Btrigger = 0, Blim = Bemp, Bpa = Bpa,
  Fscan = Fscan,
  verbose = FALSE,
  # Nrun = c(fy-iy), 
  keep.sims = TRUE
)



# subset(FIT2$rbp, variable == "Recruitment" & Ftarget == 0)$p50
# df2 |> filter(year > 2040, slot == "rec", Ftarg == 0) |> summarise(median(data)) 

### Btrigger, Flim - noErr, noAR ----

rbp <- tibble(FIT2$rbp)
Fs <- rbp |> filter(variable == "Spawning stock biomass") |> select(Ftarget)
Fs <- Fs$Ftarget
Bq05 <- rbp |> filter(variable == "Spawning stock biomass") |> select(p05)
Bq05 <- Bq05$p05
Bq50 <- rbp |> filter(variable == "Spawning stock biomass") |> select(p50)
Bq50 <- Bq50$p50

head(rbp %>% filter(variable == "Spawning stock biomass", Ftarget == 0))

plot(Fs, Bq05, ylab="tonnes", xlab="F", main = '5% percentile of SSB versus F')
abline(v = Fmsy, lty = 1, col = 2)
spl <- smooth.spline(x = Fs, y = Bq05)
lines(spl)
(Btrig <- round(predict(spl, x = Fmsy)$y))


### Bpa - based on Blim. Trumps Btrigger, if larger ----

Btrig_used <- max(Btrig, Bpa) # used in AR in 4.3


## 4.3. wErr_wAR ------------------------

set.seed(1111)
FIT3 <- eqsim_run(srr_fit, 
  bio.years = range(bio_years), 
  sel.years = range(sel_years),
  bio.const = bio_const, 
  sel.const = sel_const,
  rhologRec = mean(srpars$rho),
  Fcv = Fcv, Fphi = Fphi, SSBcv = SSBcv,
  Btrigger = Btrig_used, Blim = Bemp, Bpa = Bpa, # *** Btrig updated
  Fscan = Fscan,
  verbose = FALSE,
  # Nrun = c(fy-iy),
  keep.sims = TRUE
)

(Fpa <- round(FIT3$Refs2["catF","F05"], 3))



# eqsim_plot(FIT3)
eqsim_plot_range(FIT1)
eqsim_plot_range(FIT2)
eqsim_plot_range(FIT3)


## 4.4. summary BRPs ------------------------------------------------------------
eqsimBRPs <- BRPs*NA
eqsimBRPs[, "Fpa"] <- Fpa
eqsimBRPs[, "Fmsy"] <- Fmsy
eqsimBRPs[, "FmsyL"] <- FmsyL
eqsimBRPs[, "FmsyH"] <- FmsyU
eqsimBRPs[, "Btrigger"] <- Btrig
eqsimBRPs[, "Blim"] <- round(Bemp)
eqsimBRPs[, "Bpa"] <- Bpa
eqsimBRPs[, "MSY"] <- MSY
eqsimBRPs[, "Bmsy"] <- Bmsy
eqsimBRPs[, "B0"] <- B0
eqsimBRPs


eqsim_default <- list(
  wErr_noAR = FIT1,
  noErr_noAR = FIT2,
  wErr_wAR = FIT3)



# 5. save results ----------------------------------------------------------

save(eqsimBRPs, file = "model/eqsimBRPs.Rdata")

save(eqsim_default, eqsimBRPs, file = "model/eqsim_default_eqsimBRPs.Rdata")

