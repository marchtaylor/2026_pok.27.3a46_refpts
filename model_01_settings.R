# remove all existing objects ----------------
rm(list = ls())



# get objects -------------------------------------------------------------
load("data/objects.Rdata")

DIMS <- dims(stock)
RANGE <- range(stock)


# settings ----------------------------------------------------------------

# NUMBER of cores
ncores <- 10 # detectCores()/2 

# ITERATIONS
it <- 300

# PROJECTION years
ny <- 30

# EVALUATION years - number of terminal projections years to evaluate for reference points
ey <- 10
if(ey > ny){stop("Number of evaluation years, ey, must be smaller than numver of projection years, ny")}

# DATA year - last one
dy <- RANGE["maxyear"]

# INTERMEDIATE year
iy <- RANGE["maxyear"] + 1

# FINAL year
fy <- dy + ny

# Fscan levels
Fscan <- seq(0, 0.8, by = 0.02)

# lags
data_lag = 1
management_lag = 1

data.frame(ncores, it, ny, ey, dy, iy, fy, data_lag, management_lag)

# assessment error
Fcv <- 0.212
Fphi <- 0.423
SSBcv <- 0

# stock recruitment settings
recrLag <- RANGE["min"]
srr_years <- (2004 - recrLag) : dy

srr_models <- "segreg3"
bio_const <- FALSE
bio_years <- (dy-10):dy
sel_const <- FALSE
sel_years <- (dy-5):dy
use_rhologRec <- TRUE



save(ncores, it, ny, ey, dy, iy, fy, Fscan, 
  data_lag, management_lag,
  DIMS, RANGE, 
  Fcv, Fphi, SSBcv,
  srr_years, srr_models, bio_const, bio_years, sel_const, sel_years, use_rhologRec,
  file = "model/settings.Rdata")

