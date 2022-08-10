%macro CXO_wt(data, exposure, event, Id, B=1000);
*calculation of weights (w0 and w1) for binary exposure data;
**case-crossover study only with no time controls;
**data is the dataset; 
**exposure is the binary exposure variable;
**event is the outcome (=1 for case-crossover study), Id is the patient Id;
**B is the number of bootstrapped replicates, default is 1000;


data d17;
	set &data.;
	
	e=&exposure.;
	PtID = &Id.;
	case=&Event.;	
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
	
	dummy=1;
run;
	

*** bootstrapped standard errors;

Proc surveyselect data=pt out=pt_boots
	Seed=4321
	Method=urs
	Samprate=1
	outhits
	Rep=&B.;
run;
*pt_bs has bootstrapped replicates;

proc sort data=pt_boots out=pt_bs nodupkey; 
	by Replicate PtID NumberHits; 
run;

*initialise bootstrapped OR estimates;
data est_boot;
run;

**calculate weighted OR for every bootstrapped replicate;
%for i = 1 to &B.;
	data pt_bs;
		set pt_boots(where=(replicate EQ i));
		
	run;
	**create new ID to account for patients who are sampled mopre than once;
	data id(keep = PtID id;
		set pt_bs(where=(case Eq 1));
		id = _N_;
	run;
	
	data pt_bs;
		merge pt_bs(in=a) id(in=b);
		by PtID;
	run;

	%CXO_wt(pt_bs, exposure=e, event=case, Id=id);
	
	data est_CXO;
		set est_CXO(keep=Variable OR_G rename=(OR_G = OR_G_bs);
		
		rep=&i.;
	run;
	
	proc append base=est_boot data=est_CXO force;
	run;

%end;

*calculate bootstrapped estimate from median and CI from 2.5 and 97.5 percentiles;
proc sort data=est_boot;
	by Variable;
run;

proc univariate data=est_boot noprint; 
	var OR_G_bs; 
	by Variable;
	output out=est_CXO_boot n=n nobs=nobs pctlpts=2.5, 50, 97.5 pctlpre=OR_G;
run;
**est_CXO_boot contains the estimates for weighted CL;


%mend CXO_wt;
