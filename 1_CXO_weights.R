### create weights
rm(list=ls())

library(dplyr)
library(survival)

load(file="CXO.Rdata")
source("CXO_funcs.R")


cfit_wt <- CXO_wt(cases, exposure = ex, event = Event, Id=Id)
summary(cfit_wt)
exp(cbind(coef(cfit_wt), confint(cfit_wt)))  ##need to use bootstrap for CIs

### create weights
rm(list=ls())

library(dplyr)
library(survival)

load(file="CXO.Rdata")
source("CXO_funcs.R")


cfit_wt <- CXO_wt(cases, exposure = ex, event = Event, Id=Id)
summary(cfit_wt)
exp(cbind(coef(cfit_wt), confint(cfit_wt)))  ##need to use bootstrap for CIs

