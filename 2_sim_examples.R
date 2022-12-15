##simulated data examples
##case-time-controls
##examples 5-8 with time-varying confounder z


rm(list=ls())

library(dplyr)
library(survival)

load(file="Data/ctcsim.Rdata")
source("CXO_funcs.R")