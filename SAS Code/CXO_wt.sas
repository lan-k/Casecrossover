%macro CXO_wt(data, exposure, event, Id);
*calculation of weights (w0 and w1) for binary exposure data;
**case-crossover study only with no time controls;
**data is the dataset; expsoure is the binary exposure variable;
**event is the outcome (=1 for case-crossover study), Id is the patient Id;


data d17;
	set &data.;
	
	e=&exposure.;
	PtID = &Id.;
	case=&Event.;	
run;

data d23; 
	set d17; 
	by PtID; 
	retain caseV PT01 PT10 PT0CXO PT1CXO;
	if first.PtID then do; 
		caseV=e; 
		PT01=0; 
		PT10=0;  
		PT0CXO=ifn(e=0,1,0); 
		PT1CXO=ifn(e=1,1,0);
	end;
	else do; 
		PT01=PT01+ifn(casev=0 and E=1,1,0); 
		PT10=PT10+ifn(casev=1 and E=0,1,0);
		PT0CXO=PT0CXO+ifn(E=0,1,0); PT1CXO=PT1CXO+ifn(E=1,1,0);
	end;
run;
data d23; 
	set d23; 
	by PtID;
	if last.PtID then output;
run;
data d23; 
	set d23;
	concordant=IFN(PT0CXO=0, 0, 1)+IFN(PT1CXO=0,0,1);
run;

*exclude concordant cases;
data d23(drop=concordant); 
	set d23; 
	if concordant^=1;
run;

data d24; 
	set d23; 
	retain a0 0 a1 0 PT10m 0 PT01m 0;
	PT10m=PT10m+PT10;
	PT01m=PT01m+PT01;
	a1=a1+caseV;
	a0=a0+(1-caseV);
run;
data d24; 
	set d24 end=final;
	if final then output;
run;

data d24(keep=dummy n00 n10); 
	set d24;
	dummy=1; 
	PT10m=PT10m/a1;
	PT01m=PT01m/a0;
	*estimation of pik0: pik0 is defined as pik/pi0;
	pi00=1; 
	pi10=PT01m/PT10m;
run;

data d23; 
	set d23; 
	dummy=1;
run;

data d23; 
	merge d23 d24; 
	by dummy;
run;

data d23(drop=dummy); 
	set d23; 
	w0=pi00/PT0CXO; 
	w1=pi10/PT1CXO; 
run;

data d25(keep=PtID w1 w0); 
	set d23;
run;

proc sort data=d25; 
	by PtID;
run;

data d25; 
	merge d25(in=ina) d17(in=inb); 
	by PtID; 
	if ina=1 and inb=1;
run;

data d25; 
	set d25; 
	if e=1 then do; 
		wt=w1; 
	end;
	else do; 
		wt=w0; 
	end;
	lw=log(wt);
run;

*weighted conditional logistic regression;
proc logistic data=d25 descending; 
	model case=e /offset=lw; 
	strata PtID; 
	ods output parameterestimates=d61 oddsratios=d62;
run;


**dataset est_CXO contains the weighted odds ratios;
data est_CXO(keep=variable CL_est CL_SE CL_OR CL_OR_L CL_OR_U); 
	merge d61(in=ina keep=variable estimate stderr) d62 (in=inb rename=(Effect=Variable)); 
	by Variable; 
	if ina=1;

	CL_OR=OddsRatioEst; 
	CL_OR_L=LowerCL; 
	CL_OR_U=UpperCL;
	CL_est=Estimate; 
	CL_SE=StdErr;
run;	


%mend CXO_wt;
