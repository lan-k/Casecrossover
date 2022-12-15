libname sim "/home/lankelly0/CXO/Data";

%include "/home/lankelly0/CXO/SAS Code/CXO_wt.sas";
%include "/home/lankelly0/CXO/SAS Code/CXO_tc_wt.sas";


proc import datafile="/home/lankelly0/CXO/Data/sampdata_scenario01.csv" out=dat
		dbms=csv replace;

run;

data dat;
	set dat;
	
	event = case_period * (1-tc);
run;


data cases;
	set dat(where=(tc NE 1));
run;

%CXO_wt(cases, exposure=ex, event=event, Id = Pt_ID)
	
proc print data=out; run;


%macro test_sim(i);

	proc import datafile="/home/lankelly0/CXO/Data/sampdata_scenario0&i..csv" out=dat
		dbms=csv replace;
	run;
	
	data dat;
		set dat;
	
		event = case_period * (1-tc);
	run;
	
	%CXO_tc_wt(dat, exposure=ex, event=event, Id = Pt_ID)
	
	%put &i.;
	proc print data=out; run;


%mend;

%test_sim(1);
%test_sim(2);
%test_sim(3);
%test_sim(4);

