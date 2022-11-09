### create weights
rm(list=ls())

library(dplyr)
library(survival)

load(file="CXO.Rdata")
source("CXO_funcs.R")

#case-crossover
cfit_wt <- CXO_wt(cases, exposure = ex, event = Event, Id=Id)
summary(cfit_wt)
exp(cbind(coef(cfit_wt), confint(cfit_wt)))  ##need to use bootstrap for CIs


##with bootstrapped SEs

cfit_wt_boot <- CXO_wt_boot(cases, exposure = ex, event = Event, Id=Id, B=500)

cfit_wt_boot

#case-time-control

tcfit_wt <- CXO_tc_wt(casetimecontrols, exposure = ex, event = Event, Id=Id)
exp(cbind(coef(tcfit_wt), confint(tcfit_wt)))

(tcfit_wt_boot <- CXO_tc_wt_boot(casetimecontrols, exposure = ex, event = Event, Id=Id, B = 500))

