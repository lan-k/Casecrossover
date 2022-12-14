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


##with bootstrapped CIs
(cfit_wt_boot <- CXO_wt_boot(cases, exposure = ex, event = Event, Id=Id, B=50))

#case-time-control

tcfit_wt <- CXO_tc_wt(casetimecontrols, exposure = ex, event = Event, Id=Id)
exp(cbind(coef(tcfit_wt), confint(tcfit_wt)))
##with bootstrapped CIs
(tcfit_wt_boot <- CXO_tc_wt_boot(casetimecontrols, exposure = ex, event = Event, Id=Id, B = 50))


# with time varying confounder
# case-crossover with time varying confounder z
load(file="cases_tvc.rds")
CXO_wt(data=cases_tvc, exposure = ex, event = event, Id=Pt_ID, tvc = z)
##bootstrapped CIs
(cfit_wt_z_boot <- CXO_wt_boot(data=cases_tvc, exposure = ex, 
                                   event = event, Id=Pt_ID, tvc = z, B=50))

# case-time-control with time varying confounder z
load(file="casetimecontrols_tvc.rds")
tcfit_wt_z <- CXO_tc_wt(data=casetimecontrols_tvc, exposure = ex, event = event, Id=Pt_ID, tvc = z)
exp(cbind(coef(tcfit_wt_z), confint(tcfit_wt_z)))
##bootstrapped CIs
(tcfit_wt_z_boot <- CXO_tc_wt_boot(data=casetimecontrols_tvc, exposure = ex, 
                                  event = event, Id=Pt_ID, tvc = z, B=50))

