
%include "CXO_wt.sas";
%include "CXO_tc_wt.sas";


proc import datafile="../Data/sampdata_scenario01.csv" out=dat
		dbms=csv replace;

run;

*cases only;
data caseids;
	set dat(where=(event EQ 1));
run;

data cases;
	merge dat(in=a) caseids(in=b);
	by Pt_ID;
	if b;
run;

%CXO_wt(cases, exposure=ex, event=event, Id = Pt_ID)
	
proc print data=out; run;


*case-time-controls;
%macro test_sim(i);

	proc import datafile="../Data/sampdata_scenario0&i..csv" out=dat
		dbms=csv replace;
	run;
	

	%CXO_tc_wt(dat, exposure=ex, event=event, Id = Pt_ID)
	
	%put &i.;
	proc print data=out; run;


%mend;

%test_sim(1);
%test_sim(2);
%test_sim(3);
%test_sim(4);

