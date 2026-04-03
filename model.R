# Run analysis, write model results

# Before:
# After:

library(icesTAF)

mkdir("model")

source("model_01_settings.R")
source("model_02_propagate_stock.R")
source("model_03_srr_fits.R")
source("model_04_condition_mse.R")
source("model_05_run_mse.R")

