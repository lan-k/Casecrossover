rm(list=ls())
library(dplyr)
library(survival)
library(pubh)

load(file="drugdata.rds")  #drugdata from the WCE package

drugdata <- drugdata %>%
  group_by(Id) %>%
  mutate(futime=max(Stop),
         ex=as.numeric(dose > 0)) %>%
  ungroup()


caseids <- drugdata %>% 
  filter(Event == 1) 

caseids <- caseids %>%
  select(Id) %>%
  left_join(drugdata, by="Id")


##create a CXO study with 90 day time window

fu = 90

timecontrols <- drugdata %>%
  filter(futime >= fu)

cases <- caseids %>%
  filter(futime >= fu) %>%
  group_by(Id) %>%
  mutate(minex=min(ex), maxex=max(ex),
         concordant = minex==maxex,
         day=Stop-futime+fu,
         wt=1) %>%
  filter(!concordant, day>0) %>%
  ungroup()


save(drugdata, cases, file="CXO.Rdata")
##conditional logistic regression

cfit <- clogit(Event ~ ex + strata(Id) , data=cases) #+ offset(wt)
summary(cfit)
exp(cbind(coef(cfit), confint(cfit)))

## M-H estimate

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
  
  names(res) <- c("Common OR", "Lower CI", "Upper CI", "Pr(>|z|)")
  rownames(res) <- c()
  
  return(res)
}



mh <- mhor(formula = Event ~ Id/ex, data=cases) 
mh

##case-time-control study with 90 day time window
# %macro select_control();
# %do i=1 %to &n_obs;
# *dcase_i is one case selected from cases in turn;
# data dcase_i(keep=ite case_id day_case); set dcase; if _n_=&i; case_id=Pt_id;
# *dcontrol_c_i includes candidates of controls for each case who are at risk, and have not yet
# been selected as a control (control=0);
# data dcontrol_c_i; merge dcontrol_candidates(in=ina) dcase_i(in=inb); by ite; if ina=1 and inb=1;
# data dcontrol_c_i; set dcontrol_c_i; if day_case>=day_first and day_case<=day_last and control=0;
# *one control is selected for each case;
# data dcontrol_c_i; set dcontrol_c_i; 
# %let kk=%eval(6543210+&ite*1000+&i);
# call streaminit(&kk);
# x=rand("uniform");
# proc sort data= dcontrol_c_i; by x;
# data dcontrol_i; set dcontrol_c_i; if _n_=1;
# data dcontrol_i(keep=ite Pt_ID day_first day_case case_id); set dcontrol_i;
# *dcontrol is the list of controls for the same number of cases;
# %if &i=1 %then %do; data dcontrol; set dcontrol_i; %end;
# %else %do; data dcontrol; set dcontrol dcontrol_i; %end;
# *the value 1 is given to 'control' when selected as a control; 
# data dcontrol_candidates; merge dcontrol_candidates (in=ina) dcontrol_i(in=inb); by ite Pt_id; if ina=1;
# *control=1 is given once selected as a control;
# data dcontrol_candidates; set dcontrol_candidates; if case_id^=. then control=1;
# data dcontrol_candidates(drop=day_case case_id); set dcontrol_candidates;
# run;
# %end;
# data dcase; set dcase; case=1;
# data dcontrol; set dcontrol; case=0;
# run;
# %mend select_control;
