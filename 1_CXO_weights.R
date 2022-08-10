### create weights
rm(list=ls())

library(dplyr)
library(survival)

load(file="CXO.Rdata")
source("CXO_funcs.R")


cfit_wt <- CXO_wt(cases, exposure = ex, event = Event, Id=Id)
summary(cfit_wt)
exp(cbind(coef(cfit_wt), confint(cfit_wt)))  ##need to use bootstrap for CIs


cfit_wt_boot <- CXO_wt_boot(cases, exposure = ex, event = Event, Id=Id, B=200)

cfit_wt_boot