**weighted OR for case-time control studies;

%macro CXO_tc_wt(data, exposure, event, Id,   out=out);
** calculation of weights (w0 and w1) for binary exposure data;
** data is the dataset, assumed to be in long format with one row per Id per period; 
** exposure is the binary exposure variable ;
** event is the outcome (=1 for cases at case period, 0 for cases in control periods and 0 for time controls in all periods);
** Id is the patient Id;


**NOTE: data is assumed to sorted so that the case period is the last row per ID;

data pt;
	set &data.;
	
	e=&exposure.;
	PtID = &Id.;
	event=&Event.;	
run;


**remove discordant cases;
proc summary data=pt nway;
	by PtID;
	var e event;
	output out=discordant(drop=_FREQ_ _TYPE_) max(e)=max min(e)=min max(event) = d;
run;
*d=1 for all periods of cases, 0 for all periods of time-controls;


data discordant(where=(discordant Eq 1) drop=min max);
	set discordant;
	
	discordant = min NE max;
run;


data pt(drop=discordant);
	merge pt(in=a) discordant(in=b);
	by PtID;
	if b;
run;

data pt;
	set pt;
	by PtID;
	unex = 1-e; *unexposed;
	if last.PtID then case_period = 1;  *case period is the last period per ID for both cases and time-controls;
	else case_period = 0; 
	control_period = 1-case_period;
	c1= e EQ 1 and case_period EQ 1; *exposed case period for both cases and time controls;
	control_period_ex = e* control_period;  *exposed control period;
	control_period_unex = unex* control_period;  *unexposed control period;
run;



**calculate the weights;

proc summary data=pt nway;
	class PtID;
	types PtID;
	var c1 e unex control_period_ex control_period_unex;
	output out = dpt(drop=_TYPE_ _FREQ_) max(c1) = c1 
		sum(e unex control_period_ex control_period_unex) = PT1CXO PT0CXO control_period_ex control_period_unex;
run;

data dpt(keep = PtID c0 c1 PT01 PT10 PT1CXO PT0CXO dummy);
	set dpt;
	
	c0=1-c1;
	if c1 EQ 1 then do;  *exposed case period;
		PT10 = control_period_unex;	*sum of unexposed control periods;
		PT01 = 0;
	end;
	else do; *unexposed case period;
		PT01 = control_period_ex; *sum of exposed control periods;
		PT10 =0;	
	end;
	dummy=1;
	
run;



proc summary data=dpt nway;
	var c0 c1 PT10 PT01;
	output out = n_tc(drop=_TYPE_ _FREQ_) sum(c0 c1 PT10 PT01) = n0 n1 PT10 PT01;
run;

**n0=a0+b0 = number of cases and time controls with an unexposed case period;
**n1=a1+b1 = number of cases and time controls with an exposed case period;

data n_tc;
	set n_tc;
	
	PT01m = PT01/n0;
    PT10m = PT10/n1;
    pi00=1;
    pi10=PT01m/PT10m;
	dummy=1;
	
run;


data dpt(drop=dummy);
	merge dpt(in=a keep=PtID c0 c1 PT1CXO PT0CXO dummy) n_tc(in=b keep=pi00 pi10 dummy);
	if a;
	by dummy;
	
    w0=pi00/PT0CXO;
    w1=pi10/PT1CXO;
run;

data cases_wt;
	merge pt(in=a) dpt(in=b keep=PtID c0 c1 w0 w1);
	by PtID;
	if a;
	
	wt = ifn(e EQ 1, w1, w0);
	lw=log(wt);
	ex_tc=e;  *KK coefficient of ex_tc is for time-controls;
	ex=d*ex_tc;  *coefficient of e is for exposure-outcome association;
	
run;


*weighted conditional logistic regression for each bootstrapped sample;
proc logistic data=cases_wt descending; 
	model case_period=ex ex_tc /offset=lw; 
	strata PtID; 
	ods output  oddsratios=OR;
run;



**dataset &out. contains the weighted OR and 95% CI without bootstrapping;
data &out.(keep=variable OR_G OR_G_L OR_G_U); 
	set OR(rename=(Effect=Variable)); 

	OR_G=OddsRatioEst; 
	OR_G_L=LowerCL; 
	OR_G_U=UpperCL;
run;	



%mend CXO_tc_wt;
