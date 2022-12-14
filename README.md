# Casecrossover
This repository contains software for weighted case-crossover analysis described in 

Kubota K, Kelly TL, Sato T, Pratt N, Roughead E, Yamaguchi T. A novel weighting method to remove bias from within-subject exposure dependency in case-crossover studies. 
BMC Med Res Methodol  2021; 21:214

https://bmcmedresmethodol.biomedcentral.com/articles/10.1186/s12874-021-01408-5

Upcoming article describing methods to reduce bias due to chronic exposure with and without switching:

Kubota K, Kelly TL. Bias due to within-subject exposure dependency with or without bias due to lack of pairwise exchangeability when exposure is chronic in case-crossover and case-time-control studies: A simulation study. Am J. Epidem., in press.

**The software is currently undergoing testing, please use with caution.**

# Code
There are R and SAS versions of the code. 

R code examples: 1_CXO_weights.R
SAS macro: SAS Code/CXO_wt.sas

There are separate R functions for case-crossover (CXO_wt, CXO_wt_boot) and case-time-control (CXO_tc_wt, CXO_tc_wt_boot) designs. They are available with or without a binary time-varying confounder.

Other files are data set up scripts or functions for the weighted analysis.


# Data

The data should be structured with one row per period per person. The data should be ordered so that for each person, the earliest period appears first in the data, the case period appears last.

There should be the same number of control periods and therefore rows of data for each person.

## Inputs

The following 3 variables are required:
Patient ID
Binary exposure indicator
Binary indicator for the outcome, 0 in control periods, 1 in case period for cases

Optional:
Binary time-varying confounder

## Outputs

CXO_wt, CXO_tc_wt output a 'clogit' object.
CXO_wt_boot, CXO_tc_wt_boot output the following columns:

- Variable: Covariates in model
- est0: OR from weighted conditional logistic regression
- est: bootstrapped OR using mean if normal approximation is used, median otherwise
- lower: bootstrapped lower 95% CI from normal approximation if used, 2.5th quantile otherwise
- upper: bootstrapped upper 95% CI from normal approximation if used, 97.5th quantile otherwise

Output variables (rows) are:

- ex: exposure of cases
- ex_tc: exposure of time controls (if used)
- z: time-varying confounder of cases (if used)
- z_tc: time-varying confounder of time controls (if used)

More examples and an R package to come! Watch this space.
