### create weights
rm(list=ls())

library(dplyr)
library(survival)

load(file="CXO.Rdata")
dperiods <- cases %>% 
  group_by(Id) %>%
  mutate(case_period = Event == 1,
         control_period = Event !=1, 
         c1=as.numeric(ex==1 & case_period) #exposed case period
         ) %>%  
  summarise(c1=max(c1),
            c0=1-c1, #unexposed case period
            PT10 = ifelse(c1==1, sum((1-ex)*control_period),0),  #number of unexposed control periods per person
            PT01 = ifelse(c1==0, sum(ex*control_period),0),  #number of unexposed control periods per persons
            # 
         ) %>%  
  ungroup()

dsumm <- dperiods %>%
  mutate(a1=sum(c1), #number of exposed case periods
         a0=sum(c0))  # numer of unexposed case periods


# *PT10 is the number of unexposed control periods of the exposed case;
# *PT01 is the number of exposed control period of the unexposed case;
# PT10=IFN(c0=1 and c1=0, 1, 0)+IFN(c0=1 and c2=0, 1, 0);
# PT01=IFN(c0=0 and c1=1, 1, 0)+IFN(c0=0 and c2=1, 1, 0);
# *PT0CXO is the number of unexposed (case or control) periods;
# *PT1CXO is the number of exposed (case or control) periods;
# PT0CXO=(1-c0)+(1-c1)+(1-c2);
# PT1CXO=c0+c1+c2;
# dummy=1;
# data dpai; set dcase; retain a0 0 a1 0 PT10m 0 PT01m 0;
# PT01m=PT01m+PT01;
# PT10m=PT10m+PT10;
# a1=a1+c0;
# a0=a0+(1-c0);
# data dperiods; merge dcase(in=ina) dpai(in=inb); by dummy; if ina=1 and inb=1;
# data dperiods(drop=dummy); set dperiods; 
# *1 record is 1 period;
# *w0 and w1 are the weight per 1 unexposed and 1 exposed period, respectively, for the Greenland likelihood;
# w0=pai00/PT0CXO; 
# w1=pai10/PT1CXO; 
# *ww0 and ww1 are the weight per 1 unexposed and 1 exposed period, respectively, for the Vines and Farrington likelihood;
# ww0=ww0/PT0CXO;
# ww1=ww1/PT1CXO;
# data dperiods; set dperiods;
# case=1; ex=IFN(c0=1,1,0); wt=IFN(c0=1, w1, w0);wtv=IFN(c0=1, ww1, ww0); output;
# case=0; ex=IFN(c1=1,1,0); wt=IFN(c1=1, w1, w0);wtv=IFN(c1=1, ww1, ww0); output;
# case=0; ex=IFN(c2=1,1,0); wt=IFN(c2=1, w1, w0);wtv=IFN(c2=1, ww1, ww0); output;
# data dperiods; set dperiods; lw=log(wt); lwv=log(wtv);
# *d1 the information for the population;
# proc logistic data=dperiods descending ; model case=ex / offset=lw; strata id; 
# ods output parameterestimates=d41 oddsratios=d42;


###SAS code##
# %macro MarkovCXO(Numberv, r11v, r10v);
# data dpopulation; Number=&Numberv; r11=&r11v; r10=&r10v;
# r01=1-r11; r00=1-r10;
# pai0=r01/(r01+r10); pai1=r10/(r01+r10);
# *pxyz is the probability that exposure status is x, y and z at the case period, 1st and 2nd controls periods, respectively.;
# p100=pai0*r00*r10; p010=pai0*r10*r01; p001=pai1*r01*r00;
# *p1e1 is the probability that the case period is exposed, and N of exposed periods is 1;
# *p0e1 is the probabilities that the case period is unexposed, and N of exposed periods is 1;
# p1e1=p100;p0e1=p010+p001;
# p110=pai0*r10*r11; p101=pai1*r01*r10; p011=pai1*r11*r01;
# *p1e2 is the probability that the case period is exposed, and N of exposed periods is 2;
# *p0e2 is the probability that the case period is unexposed, and N of exposed periods is 2;
# p1e2=p110+p101; p0e2=p011;
# Npopulation=1000*(r01+r10);
# *N1e1, N0e1, N1e2, N0e2 are the expected number of cases with p1e1, p0e1, p1e2, and p0e2, respectively;
# N1e1=Npopulation*p1e1*4; N0e1=Npopulation*p0e1;
# N1e2=Npopulation*p1e2*4; N0e2=Npopulation*p0e2;
# data dcaseinfo; set dpopulation;
# *ww0 and ww1 are total weights of Vines and Farrington for the unexposed and exposed cases, respectively.;
# c0=1; c1=0; c2=0; ww1=p1e1; ww0=p0e1; N=N1e1; output;
# c0=0; c1=0; c2=1; ww1=p1e1; ww0=p0e1; N=N0e1; output;
# c0=1; c1=1; c2=0; ww1=p1e2; ww0=p0e2; N=N1e2; output;
# c0=0; c1=1; c2=1; ww1=p1e2; ww0=p0e2; N=N0e2; output;
# data dcaseinfo(keep= c0 c1 c2 ww1 ww0 N); set dcaseinfo;
# data dcaseinfo; set dcaseinfo; retain nl;
# if _n_=1 then  nf=1; else nf=nl+1; nl=nf+n-1;
# data dcase(drop=n m nf nl); set dcaseinfo;
# do m=nf to nl; id=m; output; end;
# data dcase; set dcase;
# *PT10 is the number of unexposed control periods of the exposed case;
# *PT01 is the number of exposed control period of the unexposed case;
# PT10=IFN(c0=1 and c1=0, 1, 0)+IFN(c0=1 and c2=0, 1, 0);
# PT01=IFN(c0=0 and c1=1, 1, 0)+IFN(c0=0 and c2=1, 1, 0);
# *PT0CXO is the number of unexposed (case or control) periods;
# *PT1CXO is the number of exposed (case or control) periods;
# PT0CXO=(1-c0)+(1-c1)+(1-c2);
# PT1CXO=c0+c1+c2;
# dummy=1;
# data dpai; set dcase; retain a0 0 a1 0 PT10m 0 PT01m 0;
# PT01m=PT01m+PT01;
# PT10m=PT10m+PT10;
# a1=a1+c0;
# a0=a0+(1-c0);
# data dpai; set dpai end=final;
# if final then output;
# data dpai(keep=dummy pai00 pai10); set dpai;
# dummy=1; 
# *pai00 is pai0/pai0 and pai10=pai1/pai0;
# PT10m=PT10m/a1;
# PT01m=PT01m/a0;
# pai00=1; pai10=PT01m/PT10m;
# data dperiods; merge dcase(in=ina) dpai(in=inb); by dummy; if ina=1 and inb=1;
# data dperiods(drop=dummy); set dperiods; 
# *1 record is 1 period;
# *w0 and w1 are the weight per 1 unexposed and 1 exposed period, respectively, for the Greenland likelihood;
# w0=pai00/PT0CXO; 
# w1=pai10/PT1CXO; 
# *ww0 and ww1 are the weight per 1 unexposed and 1 exposed period, respectively, for the Vines and Farrington likelihood;
# ww0=ww0/PT0CXO;
# ww1=ww1/PT1CXO;
# data dperiods; set dperiods;
# case=1; ex=IFN(c0=1,1,0); wt=IFN(c0=1, w1, w0);wtv=IFN(c0=1, ww1, ww0); output;
# case=0; ex=IFN(c1=1,1,0); wt=IFN(c1=1, w1, w0);wtv=IFN(c1=1, ww1, ww0); output;
# case=0; ex=IFN(c2=1,1,0); wt=IFN(c2=1, w1, w0);wtv=IFN(c2=1, ww1, ww0); output;
# data dperiods; set dperiods; lw=log(wt); lwv=log(wtv);
# *d1 the information for the population;
# data d1(keep=Number r11 r10 Ncases); set dpopulation ;Ncases=N1e1+N0e1+N1e2+N0e2;
# *OR_SCL standard CL;
# proc logistic data=dperiods descending; model case=ex ; strata id; 
# ods output parameterestimates=d21 oddsratios=d22;
# data d2; merge d21(in=ina) d22 (in=inb); 
# data d2(keep= OR_SCL OR_SCL_L OR_SCL_U); set d2;
# OR_SCL=OddsRatioEst; OR_SCL_L=LowerCL; OR_SCL_U=UpperCL;
# *OR_VF Vines and Farrington;
# proc logistic data=dperiods descending ; model case=ex / offset=lwv; strata id; 
# ods output parameterestimates=d31 oddsratios=d32;
# data d32(rename=(Effect=Variable)); set d32;
# data d3; merge d31(in=ina) d32 (in=inb);
# data d3; set d3; if variable="ex" then output;
# data d3(keep= OR_VF OR_VF_L OR_VF_U); set d3;
# OR_VF=OddsRatioEst; OR_VF_L=LowerCL; OR_VF_U=UpperCL;
# *OR_G Greenland;
# proc logistic data=dperiods descending ; model case=ex / offset=lw; strata id; 
# ods output parameterestimates=d41 oddsratios=d42;
# data d42(rename=(Effect=Variable)); set d42;
# data d4; merge d41(in=ina) d42 (in=inb);
# data d4; set d4; if variable="ex" then output;
# data d4(keep= OR_G OR_G_L OR_G_U); set d4;
# OR_G=OddsRatioEst; OR_G_L=LowerCL; OR_G_U=UpperCL;
# *MH method;
# proc freq data=dperiods; tables id*case*ex / cmh; 
# output mhor out=d5;
# data d5(keep= OR_MH OR_MH_L OR_MH_U); set d5; 
# OR_MH=_MHOR_; OR_MH_L=L_MHOR; OR_MH_U=U_MHOR;
# *assemble results;
# data d6; merge d1-d5; 
# %if &Numberv=1 %then %do; data CXO_results; set d6; %end;
# %else %do; data CXO_results; set CXO_results d6; %end; 
# %mend;
# *MarkovCXO(Number, r11, r10);
# %MarkovCXO(1, 0.1, 0.9);
# %MarkovCXO(2, 0.9, 0.1);
# %MarkovCXO(3, 0.1, 0.1);
# run;
