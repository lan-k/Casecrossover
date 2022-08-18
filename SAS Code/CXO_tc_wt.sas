**weighted OR for case-time control studies;

%macro CXO_tc_wt(data, exposure, event, Id, tc, out=out N, M);
*calculation of weights (w0 and w1) for binary exposure data;
**data is the dataset, assumed to be in long format with one row per Id per period; 
**exposure is the binary exposure variable ;
**event is the outcome (=1 for cases at case period, 0 for cases in control periods and for time controls in all periods), Id is the patient Id;
** the case period is assumed to to be the last period per Id ;
**in all periods tc=1 for cases and 0 for time-controls ;


data d17;
	set &data.;
	
	e=&exposure.;
	PtID = &Id.;
	case=&Event.;
	tc = &tc.;	
run;


**remove discordant cases;
proc summary data=d17 nway;
	by PtID;
	var e;
	output out=discordant(drop=_FREQ_ _TYPE_) max=max min=min;
run;

data discordant(where=(discordant Eq 1) drop=min max);
	set discordant;
	
	discordant = min NE max;
run;

proc sort data=d17; 
	by PtID descending case;
run;

data pt;
	merge d17(in=a) discordant(in=b);
	by PtID;
	if b;
	
	unex = 1-e; *unexposed;
	if last.PtID then case_period = 1;  *case period is the last period per ID for both cases and time-controls;
	else case_period = 0; 
	control_period = 1-case_period;
	c1= e EQ 1 and case EQ 1; *exposed case period;
	control_period_ex = e* control_period;  *exposed control period;
	control_period_unex = unex* control_period;  *unexposed control period;
run;



**calculate the weights;

proc summary data=pt_bs nway;
	class PtID;
	types PtID;
	var c1 e unex control_period_ex control_period_unex;
	output out = dperiods(drop=_TYPE_ _FREQ_) max(c1) = c1 
		sum(e unex control_period_ex control_period_unex) = PT1CXO PT0CXO control_period_ex control_period_unex;
run;

data dperiods(keep = PtID c0 c1 PT01 PT10 PT1CXO PT0CXO);
	set dperiods;
	
	c0=1-c1;
	if c0 EQ 1 then do;  *unexposed case period;
		PT10 = control_period_unex;	
		PT01 = 0;
	end;
	else do; *exposed case period;
		PT01 = control_period_ex;
		PT10 =0;	
	end;
run;



proc summary data=dperiods nway;
	var c0 c1 PT10 PT10;
	output out = n_tc(drop=_TYPE_ _FREQ_) sum(c0 c1 PT10 PT10) = n0 n1 PT10 PT10;
run;

**n0=a0+b0 = number of cases and time controls with an unexposed case period;
**n0=a1+b1 = number of cases and time controls with an exposed case period;

data dperiods;
	merge dperiods(in=a) n_tc(in=b);
	if a;
	
	PT01m = PT01/n0
    PT10m = PT10/n1
    pi00=1
    pi10=PT01m/PT10m
    w0=pi00/PT0CXO
    w1=pi10/PT1CXO
run;

data cases_wt;
	merge pt(in=a) dperiods(in=b keep=PtID c0 w0 w1);
	by PtID;
	if a;
	
	wt = ifn(e EQ 1, w1, w0);
	lw=log(wt);
	ex_tc=tc*e;
run;


*weighted conditional logistic regression for each bootstrapped sample;
proc logistic data=cases_wt descending; 
	model case=e ex_tc /offset=lw; 
	strata PtID; 
	ods output  oddsratios=OR(rename=(Effect=Variable));
run;



**dataset &out. contains the weighted OR and 95% CI without bootstrapping;
data &out.(keep=variable CL_est CL_SE OR_G OR_G_L OR_G_U); 
	set OR(rename=(Effect=Variable)); 

	OR_G=OddsRatioEst; 
	OR_G_L=LowerCL; 
	OR_G_U=UpperCL;
run;	



%mend CXO_tc_wt;
