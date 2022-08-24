%macro CXO_wt_boot(data, exposure, case, Id, B=1000, seed=4321, out=out);
*calculation of weights (w0 and w1) for binary exposure data;
**case-crossover study only with no time controls;
**data is the dataset; 
**exposure is the binary exposure variable;
**case is the outcome (=1 for case-crossover study), Id is the patient Id;
**B is the number of bootstrapped replicates, default is 1000;


data pt;
	set &data.;
	
	e=&exposure.;
	PtID = &Id.;
	event=&Event.;	
	
	unex = 1-e; *unexposed;
	case_period  = event EQ 1;  *assuming no time controls, case period is where the event occurs;
	control_period = 1-case_period;
	c1= e EQ 1 and event EQ 1; *exposed case period;
	control_period_ex = e* control_period;  *exposed control period;
	control_period_unex = unex* control_period;  *unexposed control period;
run;


**remove discordant cases;
proc summary data=pt nway;
	by PtID;
	var e;
	output out=discordant(drop=_FREQ_ _TYPE_) max=max min=min;
run;

data discordant(where=(discordant Eq 1) drop=min max);
	set discordant;
	
	discordant = min NE max;
run;


data pt(drop=discordant);
	merge pt(in=a) discordant(in=b);
	by PtID;
	if b;
	
run;


**create dataset of Ids for resampling;
data id(keep=PtID);
	set pt;
	by PtID;
	if last.PtID then output;
run;
	

*** bootstrapped standard errors;

Proc surveyselect data=id noprint out=pt_bs
	Seed=&seed.
	Method=urs
	Samprate=1
	outhits
	Rep=&B.;
run;
*pt_bs has bootstrapped replicates of ids;


data pt_bs;
	set pt_bs;
	newid = _N_;  * create new patient ID;
run;


data pt_bs;
		merge pt_bs(in=a keep=replicate PtID newid) pt(in=b);
		by PtID;
		if a;			
run;

**calculate the weights;

proc summary data=pt_bs nway;
	class replicate newid;
	types replicate*newid;
	var c1 e unex control_period_ex control_period_unex;
	output out = dpt(drop=_TYPE_ _FREQ_) max(c1) = c1 
		sum(e unex control_period_ex control_period_unex) = PT1CXO PT0CXO control_period_ex control_period_unex;
run;

data dpt(keep = replicate newid c0 c1 PT01 PT10 PT1CXO PT0CXO);
	set dpt;
	
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


proc summary data=dpt nway;
	class replicate;
	types replicate;
	var c0 c1 PT10 PT10;
	output out = n_case(drop=_TYPE_ _FREQ_) sum(c0 c1 PT10 PT10) = a0 a1 PT10 PT10;
run;

dta n_case;
	set n_case;
	
	PT01m = PT01/a0;
    PT10m = PT10/a1;
    pi00=1;
    pi10=PT01m/PT10m;
	
run;

data dpt;
	merge dpt(in=a) n_case(in=b);
	by replicate;
	if a;
		
    w0=pi00/PT0CXO;
    w1=pi10/PT1CXO;
run;

data cases_wt;
	merge pt_bs(in=a) dpt(in=b);
	by replicate newid;
	if a;
	
	wt = ifn(e EQ 1, w1, w0);
	lw=log(wt);
run;
	
	
*weighted conditional logistic regression for each bootstrapped sample;
proc logistic data=cases_wt descending; 
	by replicate;
	model event=e /offset=lw; 
	strata newid; 
	ods output  oddsratios=est_CXO_rep(rename=(Effect=Variable));
run;	
	
proc univariate data=est_CXO_rep(keep=replicate variable OddsRatioEst) noprint; 
	var OddsRatioEst; 
	by Variable;
	output out=est_CXO_boot n=n nobs=nobs pctlpts=50 2.5 97.5 
			pctlpre=OR_G pctlname = _median _L _U;
run;
**est_CXO_boot contains the estimates for weighted CL;

**calculate nonbootstrapped estimate;
%CXO_wt(pt, exposure=e, case=case, Id=newid, out=est_CXO_orig);


**combine the nonbootstrapped estimate with the bootstrapped;
data &out.;
	merge est_CXO_boot(in=a) est_CXO_orig(in=b keep=Variable OR_G);
	by Variable;
	if a;
run;


%mend CXO_wt_boot;
