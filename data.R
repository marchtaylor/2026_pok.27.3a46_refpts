# Prepare data, write CSV data tables

# Before:
# After:

library(icesTAF)

mkdir("data")


# load packages -----------------------------------------------------------
library(FLCore)
library(FLasher)
library(mse)
library(msemodules)
library(FLSRTMB)
library(tidyr)
library(dplyr)
library(ggplot2)
library(patchwork)
library(parallel)
library(rmarkdown)
taf.library(FLfse)
library(beepr)
library(stockassessment)


source("data_01_objects.R")

