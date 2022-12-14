###CXO functions
library(dplyr)
library(survival)
library(boot)
library(broom)
library(scales)



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
    group_by(Id, .drop = T) %>%
    rename(ex={{exposure}},
           Event={{event}},
           Id={{Id}}) %>%
    mutate(
      minex=min(ex),
      maxex=max(ex)) %>%
    filter(minex != maxex)  %>%  #remove concordant cases
    ungroup()
  
  
  cfit <- clogit(Event ~ ex + strata(Id) , data=cases, method="efron") #+ offset(wt)
  
  est_scl <- as.numeric(exp(coef(cfit)))
  
  mh <- mhor(formula = Event ~ Id/ex, data=cases)
  mh_OR <- as.numeric(mh$OR)
  
  
  bias <- percent(abs((as.numeric(est_scl) - mh_OR)/mh_OR),accuracy = 0.1)

  
  return(list( scl = est_scl, mh_or=mh_OR, bias=bias))
  
}




CXO_wt <- function(data, exposure, event, Id, tvc = NULL) {
  cases <- data %>% 
    mutate(ex={{exposure}},
           Event={{event}},
           Id={{Id}}) %>%
    group_by(Id) %>%
    mutate(
           minex=min(ex), 
           maxex=max(ex),
           case_period = Event == 1,
           control_period = Event !=1, 
           c1=as.numeric(ex==1 & case_period) #exposed case period
    ) %>%  
    filter(minex != maxex)  %>%  #remove concordant cases
    ungroup()
  
  
  dpt <- cases %>%
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
           w1=pi10/PT1CXO) %>%
    ungroup()
  
  cases_wt <- left_join(cases, dpt %>% select(Id,w0,w1), by="Id") %>% # number of unexposed case periods
    # rowwise() %>%
    ##calculate weights depending on whether period is exposed or unexposed
    mutate(wt=ifelse(ex==1, w1, w0),
           lw=log(wt),
           ex=relevel(factor(ex), ref = "0")) %>%  
    ungroup() 
  
  if (!is.null(substitute(tvc)) ) {
    cases_wt <- cases_wt %>%
      mutate(z = {{tvc}},
             z=relevel(factor(z), ref = "0")) %>%  
      ungroup() 
    wfit <- clogit(case_period ~ ex + z + strata(Id) + offset(lw) ,
                   data=cases_wt, method="efron")
  } else {
    
    wfit <- clogit(case_period ~ ex  + strata(Id) + offset(lw) , 
                  data=cases_wt, method="efron")
  }
  
  
  return(wfit)
  
}



.CI_boot <- function(data,ii, exposure, event , Id, tvc = NULL) {
  ##internal function for bootstrapping the CIs
  
  df <- data %>%
    mutate(Id= {{Id}}) 
  
  cases <- df %>% 
    group_by(Id) %>%
    slice_tail(n=1) %>%
    ungroup() %>% 
    select(Id)

  #resample the cases
  dd <- cases[ii,] %>%
    mutate(newid = row_number()) %>%
    left_join(df, by="Id")
  if (!is.null(substitute(tvc))) {
    cfit <- CXO_wt(dd, exposure = {{exposure}}, event = {{event}},tvc={{tvc}}, Id=newid)
  } else {
    cfit <- CXO_wt(dd, exposure = {{exposure}}, event = {{event}}, Id=newid)
  }  
  return(coef(cfit))
}



CXO_wt_boot <- function(data, exposure, event, Id, tvc = NULL, B=500, normal = T) {
  
  if (!is.null(substitute(tvc))) {
    df <- data %>% 
      mutate(ex={{exposure}},
             Event={{event}},
             Id={{Id}},
             z={{tvc}}) %>%
      select(ex, Event, Id, z)
    
    fitboot <- boot(data=df, statistic = .CI_boot, 
                    exposure = ex, event = Event, Id=Id, tvc=z,
                    R=B)
    
  } else {
    df <- data %>% 
      mutate(ex={{exposure}},
             Event={{event}},
             Id={{Id}}) %>%
      select(ex, Event, Id)
    
    fitboot <- boot(data=df, statistic = .CI_boot, 
                    exposure = ex, event = Event, Id=Id, 
                    R=B)
   
  }
  nvars <- dim(fitboot$t)[2]
  ci <- list()
  for (i in (1:nvars)) {
    if (normal) {
      mean = mean(fitboot$t[,i])
      sd = sd(fitboot$t[,i])
      est = exp(mean)
      lower = exp(mean - 1.96*sd)
      upper=  exp(mean + 1.96*sd)
      ci[[i]] <- data.frame(est, lower, upper)
    } else {
       temp <- quantile(fitboot$t[,i],p=c(0.5, 0.025,0.975))
       ci[[i]] <- data.frame(est=exp(temp[1]), lower=exp(temp[2]), upper=exp(temp[3]))
    }
    
    
  }
  
                   
  est0=data.frame(est0= exp(fitboot$t0))
  est0$Variable = rownames(est0)
  est <- bind_cols(est0= est0, bind_rows(ci)) %>%
    select(Variable,est0, est, lower, upper )

  return(est)
  
}


CXO_tc_wt <- function(data, exposure, event, Id, tvc = NULL) {
# tc = 1 for time-controls, 1 for cases in all periods
# data is assumed to be sorted so that the case period is the last row per Id
  case_tc <- data %>% 
    mutate(ex={{exposure}},
           Event={{event}},
           Id={{Id}}) %>%
    group_by(Id) %>%
    mutate(period = row_number(),
      minex=min(ex), 
      maxex=max(ex),
      d=max(Event), #d=1 for cases and 0 for time-controls
      case_period = as.numeric(period == max(period)),
      control_period = 1- case_period, 
      c1=as.numeric(ex==1 & case_period) #exposed case period
    ) %>%  
    filter(minex != maxex)  %>%  #remove concordant cases
    ungroup()
  
  
  dpt <- case_tc %>%
    group_by(Id) %>%
    summarise(c1=max(c1),
              c0=1-c1, #unexposed case period
              PT10 = ifelse(c1==1, sum((1-ex)*control_period),0),  #number of unexposed control periods per person
              PT01 = ifelse(c1==0, sum(ex*control_period),0),  #number of unexposed control periods per persons
              PT0CXO = sum(1-ex, na.rm=T), # number of unexposed (case or control) periods
              PT1CXO = sum(ex, na.rm=T) # number of exposed (case or control) periods
    ) %>%  
    mutate(n1=sum(c1, na.rm=T), #number of exposed case periods for both cases and time-controls
           n0=sum(c0, na.rm=T), # number of unexposed case periods for both cases and time-controls
           PT01m = sum(PT01, na.rm=T)/n0,
           PT10m = sum(PT10, na.rm=T)/n1,
           pi00=1,
           pi10=PT01m/PT10m,
           w0=pi00/PT0CXO,
           w1=pi10/PT1CXO) %>%
    ungroup()
  
  cases_wt <- left_join(case_tc, dpt %>% select(Id, w0,w1), by="Id") %>% # number of unexposed case periods
    # rowwise() %>%
    ##calculate weights depending on whether period is exposed or unexposed
    mutate(wt=ifelse(ex==1, w1, w0),
           lw=log(wt),
           ex_tc=relevel(factor(ex), ref="0"),
           ex=relevel(factor(d*ex), ref="0")) %>%  
    ungroup() 
  
  if (!is.null(substitute(tvc)) ) {
    cases_wt <- cases_wt %>%
      mutate(z_tc={{tvc}},
             z=relevel(factor(d*z_tc), ref = "0"),
             z_tc=relevel(factor(z_tc), ref = "0")) %>%  
      ungroup() 
    wfit <- clogit(case_period ~ ex  + ex_tc + z + z_tc + strata(Id) + offset(lw) , 
                   data=cases_wt, method="efron")
  } else {
    
    wfit <- clogit(case_period ~ ex + ex_tc + strata(Id) + offset(lw) , 
                   data=cases_wt, method="efron")
  }
  
  #coefficient of ex_tc is for time-controls/time-trend
  #coefficient of ex is for exposure-outcome
  
  return(wfit)
  
}

.CI_tc_boot <- function(data,ii, exposure, event , Id, tvc= NULL) {
  ##internal function for bootstrapping the CIs
  
  df <- data %>%
    mutate(Id= {{Id}}) 
  
  cases <- df %>% 
    group_by(Id) %>%
    slice_tail(n=1) %>%
    ungroup() %>% 
    select(Id)
  
  #resample the cases
  dd <- cases[ii,] %>%
    mutate(newid = row_number()) %>%
    left_join(df, by="Id")

  if (!is.null(substitute(tvc))) {
    cfit <- CXO_tc_wt(dd, exposure = {{exposure}}, event = {{event}},tvc={{tvc}}, Id=newid)
  } else {
    cfit <- CXO_tc_wt(dd, exposure = {{exposure}}, event = {{event}}, Id=newid)
  }  
                      
  
  return(coef(cfit))
}


CXO_tc_wt_boot <- function(data, exposure, event, Id, tvc=NULL, B=500, normal=T) {
  
  if (!is.null(substitute(tvc))) {
    df <- data %>% 
      mutate(ex={{exposure}},
             Event={{event}},
             Id={{Id}},
             z={{tvc}}) %>%
      select(ex, Event, Id, z)
    
    fitboot <- boot(data=df, statistic = .CI_tc_boot, 
                    exposure = ex, event = Event, Id=Id, tvc=z,
                    R=B)
  } else {
    df <- data %>% 
      mutate(ex={{exposure}},
             Event={{event}},
             Id={{Id}}) %>%
      select(ex, Event, Id)
    
    fitboot <- boot(data=df, statistic = .CI_tc_boot, 
                    exposure = ex, event = Event, Id=Id, 
                    R=B)
  }
  
  nvars <- dim(fitboot$t)[2]
  ci <- list()
  for (i in (1:nvars)) {
    if (normal) {
      mean = mean(fitboot$t[,i])
      sd = sd(fitboot$t[,i])
      est = exp(mean)
      lower = exp(mean - 1.96*sd)
      upper=  exp(mean + 1.96*sd)
      ci[[i]] <- data.frame(est, lower, upper)
    } else {
      temp <- quantile(fitboot$t[,i],p=c(0.5, 0.025,0.975))
      ci[[i]] <- data.frame(est=exp(temp[1]), lower=exp(temp[2]), upper=exp(temp[3]))
    }
  }
  
  
  est0=data.frame(est0= exp(fitboot$t0))
  est0$Variable = rownames(est0)
  est <- bind_cols(est0= est0, bind_rows(ci)) %>%
    select(Variable,est0, est, lower, upper )
  
  return(est)
  
}





