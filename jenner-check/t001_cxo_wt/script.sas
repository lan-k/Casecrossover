/* CXO_wt macro from lan-k/Casecrossover, SAS Code/CXO_wt.sas (verbatim). */
/* Weighted conditional logistic regression for weighted case-crossover analysis. */
/* Sample data below are rows drawn from the repo's Data/sampdata_scenario01.csv. */

%macro CXO_wt(data, exposure, event, Id, out=out);
*calculation of weights (w0 and w1) for binary exposure data;
**case-crossover study only with no time controls;
**data is the dataset, assumed to be in long format with one row per Id per period; 
**exposure is the binary exposure variable ;
**event is the outcome (=1 for case-crossover study), Id is the patient Id;
** the case period is assumed to to be the last period per id indicated by event=1 for CXO studies (without time controls);

data d17;
	set &data.;
	
	e=&exposure.;
	PtID = &Id.;
	event=&Event.;	
run;


proc sort data=d17; 
	by PtID descending event;
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

data d24(keep=dummy pi00 pi10); 
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
	model event=e /offset=lw; 
	strata PtID; 
	ods output parameterestimates=d61 oddsratios=d62;
run;


**dataset &out. contains the weighted odds ratios;
data &out.(keep=variable CL_est CL_SE OR_G OR_G_L OR_G_U where=(Variable NE "lw")); 
	merge d61(in=ina keep=variable estimate stderr) d62 (in=inb rename=(Effect=Variable)); 
	by Variable; 
	if ina=1;

	OR_G=OddsRatioEst; 
	OR_G_L=LowerCL; 
	OR_G_U=UpperCL;
	CL_est=Estimate; 
	CL_SE=StdErr;
run;	


%mend CXO_wt;

/* --- caller: read the sampled case-crossover data and run the macro --- */
data cases;
  input Pt_ID ex event;
  datalines;
1 1 0
1 0 0
1 0 0
1 0 0
1 0 0
1 0 0
1 0 0
1 0 0
1 0 0
1 0 0
1 0 1
2 0 0
2 0 0
2 0 0
2 0 0
2 0 0
2 0 0
2 0 0
2 1 0
2 1 0
2 1 0
2 1 1
5 0 0
5 0 0
5 1 0
5 1 0
5 1 0
5 1 0
5 1 0
5 1 0
5 1 0
5 1 0
5 1 1
6 0 0
6 0 0
6 0 0
6 1 0
6 1 0
6 1 0
6 1 0
6 1 0
6 1 0
6 1 0
6 1 1
7 0 0
7 0 0
7 0 0
7 0 0
7 0 0
7 0 0
7 0 0
7 0 0
7 0 0
7 1 0
7 1 1
11 1 0
11 0 0
11 0 0
11 0 0
11 0 0
11 0 0
11 0 0
11 0 0
11 0 0
11 0 0
11 0 1
13 0 0
13 0 0
13 0 0
13 0 0
13 0 0
13 0 0
13 0 0
13 0 0
13 0 0
13 1 0
13 1 1
15 0 0
15 0 0
15 0 0
15 1 0
15 1 0
15 1 0
15 1 0
15 1 0
15 1 0
15 1 0
15 1 1
17 1 0
17 1 0
17 0 0
17 0 0
17 0 0
17 0 0
17 0 0
17 0 0
17 0 0
17 0 0
17 0 1
19 0 0
19 0 0
19 0 0
19 1 0
19 1 0
19 1 0
19 1 0
19 1 0
19 1 0
19 1 0
19 1 1
20 0 0
20 0 0
20 0 0
20 0 0
20 0 0
20 0 0
20 0 0
20 1 0
20 1 0
20 1 0
20 1 1
23 1 0
23 1 0
23 1 0
23 0 0
23 0 0
23 0 0
23 0 0
23 0 0
23 0 0
23 0 0
23 0 1
25 1 0
25 1 0
25 1 0
25 1 0
25 1 0
25 1 0
25 1 0
25 1 0
25 1 0
25 1 0
25 0 1
26 0 0
26 0 0
26 0 0
26 0 0
26 0 0
26 0 0
26 1 0
26 1 0
26 1 0
26 1 0
26 1 1
29 0 0
29 0 0
29 0 0
29 0 0
29 1 0
29 1 0
29 1 0
29 1 0
29 1 0
29 1 0
29 1 1
31 0 0
31 0 0
31 0 0
31 0 0
31 0 0
31 1 0
31 1 0
31 1 0
31 1 0
31 1 0
31 1 1
33 0 0
33 0 0
33 0 0
33 0 0
33 0 0
33 1 0
33 1 0
33 1 0
33 1 0
33 1 0
33 1 1
34 0 0
34 0 0
34 0 0
34 0 0
34 0 0
34 0 0
34 0 0
34 0 0
34 1 0
34 1 0
34 1 1
37 1 0
37 1 0
37 1 0
37 0 0
37 0 0
37 0 0
37 0 0
37 0 0
37 0 0
37 0 0
37 0 1
38 0 0
38 0 0
38 0 0
38 1 0
38 1 0
38 1 0
38 1 0
38 1 0
38 1 0
38 1 0
38 1 1
41 0 0
41 0 0
41 0 0
41 1 0
41 1 0
41 1 0
41 1 0
41 1 0
41 1 0
41 1 0
41 1 1
43 0 0
43 0 0
43 0 0
43 0 0
43 0 0
43 0 0
43 0 0
43 0 0
43 0 0
43 1 0
43 1 1
45 1 0
45 1 0
45 1 0
45 1 0
45 1 0
45 1 0
45 1 0
45 0 0
45 0 0
45 0 0
45 0 1
46 1 0
46 1 0
46 1 0
46 0 0
46 0 0
46 0 0
46 0 0
46 0 0
46 0 0
46 0 0
46 0 1
47 1 0
47 1 0
47 1 0
47 1 0
47 1 0
47 1 0
47 0 0
47 0 0
47 0 0
47 0 0
47 0 1
51 0 0
51 0 0
51 0 0
51 0 0
51 0 0
51 0 0
51 0 0
51 1 0
51 1 0
51 1 0
51 1 1
53 0 0
53 1 0
53 1 0
53 1 0
53 1 0
53 1 0
53 1 0
53 1 0
53 1 0
53 1 0
53 1 1
54 0 0
54 0 0
54 0 0
54 0 0
54 0 0
54 0 0
54 0 0
54 0 0
54 1 0
54 1 0
54 1 1
57 0 0
57 0 0
57 0 0
57 1 0
57 1 0
57 1 0
57 1 0
57 1 0
57 1 0
57 1 0
57 1 1
59 1 0
59 1 0
59 1 0
59 0 0
59 0 0
59 0 0
59 0 0
59 0 0
59 0 0
59 0 0
59 0 1
60 0 0
60 0 0
60 0 0
60 0 0
60 0 0
60 0 0
60 0 0
60 0 0
60 1 0
60 1 0
60 1 1
63 0 0
63 0 0
63 0 0
63 0 0
63 0 0
63 0 0
63 0 0
63 0 0
63 1 0
63 1 0
63 1 1
65 0 0
65 0 0
65 0 0
65 0 0
65 0 0
65 0 0
65 1 0
65 1 0
65 1 0
65 1 0
65 1 1
67 0 0
67 0 0
67 0 0
67 0 0
67 1 0
67 1 0
67 1 0
67 1 0
67 1 0
67 1 0
67 1 1
69 0 0
69 0 0
69 0 0
69 0 0
69 0 0
69 0 0
69 0 0
69 0 0
69 0 0
69 1 0
69 1 1
71 0 0
71 0 0
71 0 0
71 0 0
71 0 0
71 0 0
71 0 0
71 0 0
71 1 0
71 1 0
71 1 1
73 0 0
73 0 0
73 0 0
73 1 0
73 1 0
73 1 0
73 1 0
73 1 0
73 1 0
73 1 0
73 1 1
74 0 0
74 0 0
74 0 0
74 1 0
74 1 0
74 1 0
74 1 0
74 1 0
74 1 0
74 1 0
74 1 1
75 0 0
75 0 0
75 0 0
75 0 0
75 0 0
75 0 0
75 1 0
75 1 0
75 1 0
75 1 0
75 1 1
76 0 0
76 0 0
76 0 0
76 0 0
76 0 0
76 0 0
76 1 0
76 1 0
76 1 0
76 1 0
76 1 1
;
run;

%CXO_wt(cases, exposure=ex, event=event, Id=Pt_ID)

proc print data=out; run;
