###CXO functions
library(dplyr)
library(survival)
library(boot)
library(broom)

CXO_wt <- function(data, exposure, event, Id) {
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



.CI_boot <- function(data,ii, exposure, event , Id) {
  ##internal function for bootstrapping the CIs
  
  cases <- data %>%
    rename(ex={{exposure}},
           Event={{event}},
           Id={{Id}}) %>%
    group_by(Id) %>%
    slice_tail(n=1) %>%
    ungroup() %>%
    select(Id)


  #resample the cases
  dd <- cases[ii,] %>%
    mutate(newid = row_number()) %>%
    left_join(data, by="Id")
  
  cfit <- CXO_wt(dd, exposure = ex, event = Event, Id=newid)
  return(coef(cfit))
}



CXO_wt_boot <- function(data, exposure, event, Id, B=500) {
  
  
  fitboot <- boot(data=data, statistic = .CI_boot, 
                  exposure = {{exposure}}, event = {{event}}, Id={{Id}},
                  R=B)
  ci <- quantile(fitboot$t,p=c(0.025,0.975))
  
  est <-data.frame(est=fitboot$t0, lower=ci[1], upper=ci[2])
  return(est)
  
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






rd <- function(data,ii){
  
  dd <- data[ii,]  #for bootstrap sample
  
  cfit<- coxph(Surv(age_admit_years, age_to_death_years_c, died_c) ~ 
                 sex  + postTH2  + stroke_type   + stroke_type*postTH2
               + estimated_prestroke_mrs   + IRSAD_quint  
               + admitted_to_stroke_unit + can_walk #+ nihss_initial + prior_stroke + prior_tia
               + ever_smoker  + pspline(n_comorbid, n=5)  #+ season
               + Hospitalname
               + strata(stroke_num) # + cluster(patid)
               , data = dd, model=T, x=T,y=T)
  
  T0 <- dd; T1 <- dd;  T2 <- dd
  
  #replace values of postTH2 
  T0$postTH2 <- 'Dec2016-Mar2017' 
  T1$postTH2 <- 'Mar2017-Sep2017'
  T2$postTH2 <- 'Sep2017-Dec2018'
  
  p0 <- predict(cfit, newdata=T0, type="expected",se.fit=TRUE,reference="strata") #expected number of events
  p1 <- predict(cfit, newdata=T1, type="expected",se.fit=TRUE, reference="strata")
  p2 <- predict(cfit, newdata=T2, type="expected",se.fit=TRUE,reference="strata")
  
  e <- dd %>% 
    mutate(expected0=p0$fit,
           expected1=p1$fit,
           expected2=p2$fit,
           surv0=1/exp(expected0),  #survival probabilities
           surv1=1/exp(expected1),
           surv2=1/exp(expected2)) 
  
  
  dr <- absrate(e) 
  
  # isch1_1_0 <- dr  %>% select(r1_0)
  # isch2_1_0 <- dr %>% filter(stroke_type == "Isch", stroke_num==2) %>% select(r1_0)
  # isch1_2_0 <- dr %>% filter(stroke_type == "Isch", stroke_num==1) %>% select(r2_0)
  # isch2_2_0 <- dr %>% filter(stroke_type == "Isch", stroke_num==2) %>% select(r2_0)
  # 
  # ICH1_1_0 <- dr %>% filter(stroke_type == "ICH", stroke_num==1) %>% select(r1_0)
  # ICH2_1_0 <- dr %>% filter(stroke_type == "ICH", stroke_num==2) %>% select(r1_0)
  # ICH1_2_0 <- dr %>% filter(stroke_type == "ICH", stroke_num==1) %>% select(r2_0)
  # ICH2_2_0 <- dr %>% filter(stroke_type == "ICH", stroke_num==2) %>% select(r2_0)
  # 
  
  # res <- c(isch1_1_0, isch2_1_0, isch1_2_0,isch2_2_0,
  #          ICH1_1_0, ICH2_1_0, ICH1_2_0, ICH2_2_0)
  
  
  
  return(c(dr$r1_0, dr$r2_0))
  
}









