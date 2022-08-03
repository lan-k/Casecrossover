3###CXO functions

CXO_wt <- function(data, exposure, event, Id) {
  cases <- data %>% 
    group_by(Id) %>%
    rename(ex={{exposure}},
           Event={{event}},
           Id={{Id}}) %>%
    mutate(
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


mhor <- function(formula, data, digits=2)  {
  #formula is Outcome ~ strata/exposure
  
  dformat <- paste('%.', digits, 'f', sep='')
  vars <- all.vars(formula)
  outcome <- vars[1]
  exposure <- vars[3]
  stratum <- vars[2]
  mht <- mantelhaen.test(data[[outcome]], data[[exposure]], 
                         data[[stratum]])
  res <- data.frame(sprintf(dformat,round(mht$estimate[[1]], digits)), 
                    sprintf(dformat,round(mht$conf.int[1], digits)), 
                    sprintf(dformat,round(mht$conf.int[2], digits)),
                    format.pval(mht$p.value, digits=2, eps=0.001))
  
  names(res) <- c("OR", "Lower CI", "Upper CI", "Pr(>|z|)")
  rownames(res) <- c()
  
  return(res)
}

SCL_bias <- function(data, exposure, event, Id) {
  
  cases <- data %>% 
    group_by(Id) %>%
    rename(ex={{exposure}},
           Event={{event}},
           Id={{Id}}) %>%
    mutate(
           minex=min(ex), 
           maxex=max(ex)) %>%  
    filter(minex != maxex)  %>%  #remove concordant cases
    ungroup()
  
  
  cfit <- clogit(Event ~ ex + strata(Id) , data=cases) #+ offset(wt)

  est_scl <- exp(coef(cfit))
  
  mh <- mhor(formula = Event ~ Id/ex, data=cases) 
  mh_OR <- as.numeric(mh$OR)
  
  
  bias <- abs(100*( est_scl - mh_OR)/mh_OR)
  
  sprintf("Bias is %.1f%%", bias)
  
  return(list(est_scl, mh_OR))
  
}


