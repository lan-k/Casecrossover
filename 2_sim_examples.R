##simulated data examples
##case-time-controls
##examples 5-8 with time-varying confounder z


rm(list=ls())

library(dplyr)
library(survival)

load(file="Data/ctcsim.Rdata")
source("CXO_funcs.R")


OR_est <- function(df, i) {
  if (i < 5) {
    temp = CXO_tc_wt(data=df, exposure = ex, event = event, Id=pt_id)
  } else {
    temp = CXO_tc_wt(data=df, exposure = ex, event = event, Id=pt_id, tvc = z)
  }
  
  return(round(exp(temp$coefficients), digits=3))
}

orlist <-list()
for (i in 1:8) {
  orlist[[i]] <- OR_est(simdata[[i]],i)
}


knitr::kable(orlist)
