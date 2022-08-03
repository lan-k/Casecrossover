### create weights
rm(list=ls())

library(dplyr)
library(survival)

load(file="CXO.Rdata")

CXO_wt <- function(data, exposure, event, Id) {
  cases <- data %>% 
    group_by(Id) %>%
    rename(ex={{exposure}},
           Event={{event}},
           Id={{Id}}) %>%
    mutate(ex={{exposure}},
           minex=min(ex), 
           maxex=max(ex),
           case_period = Event == 1,
           control_period = Event !=1, 
           c1=as.numeric(ex==1 & case_period) #exposed case period
    ) %>%  
    filter(minex != maxex)  %>%  #remove concordant cases
    ungroup()
  
  
  dperiods <- cases %>%
    group_by(Id) %>%
    summarise(c1=max(c1),
              c0=1-c1, #unexposed case period
              PT10 = ifelse(c1==1, sum((1-ex)*control_period),0),  #number of unexposed control periods per person
              PT01 = ifelse(c1==0, sum(ex*control_period),0),  #number of unexposed control periods per persons
              PT0CXO = sum(1-ex, na.rm=T), # number of unexposed (case or control) periods
              PT1CXO = sum(ex, na.rm=T) # number of exposed (case or control) periods
    ) %>%  
    mutate(a1=sum(c1, na.rm=T), #number of exposed case periods
           a0=sum(c0, na.rm=T), # number of unexposed case periods
           PT01m = sum(PT01, na.rm=T)/a0,
           PT10m = sum(PT10, na.rm=T)/a1,
           pi00=1,
           pi10=PT01m/PT10m,
           w0=pi00/PT0CXO,
           w1=pi10/PT1CXO,) %>%
    ungroup()
  
  cases_wt <- left_join(cases, dperiods %>% select(Id,w0,w1), by="Id") %>% # number of unexposed case periods
    rowwise() %>%
    ##calculate weights depending on whether period is exposed or unexposed
    mutate(wt=ifelse(ex==1, w1, w0),
           lw=log(wt)) %>%  
    ungroup() 
  
  wfit <- clogit(Event ~ ex + strata(Id) + offset(lw), data=cases_wt)
  
  return(wfit)
  
}

cfit_wt <- CXO_wt(cases, exposure = ex, event = Event, Id=Id)
summary(cfit_wt)
exp(cbind(coef(cfit_wt), confint(cfit_wt)))  ##need to use bootstrap for CIs


