# Casecrossover
Software for weighted case-crossover analysis described in 

Kubota K, Kelly TL, Sato T, Pratt N, Roughead E, Yamaguchi T. A novel weighting method to remove bias from within-subject exposure dependency in case-crossover studies. 
BMC Med Res Methodol  2021; 21:214

https://bmcmedresmethodol.biomedcentral.com/articles/10.1186/s12874-021-01408-5


There are R and SAS versions of the code. 

R code examples: 1_CXO_weights.R
SAS macro: SAS Code/CXO_wt.sas

Other files are data set up scripts or functions for the weighted analysis.


# Data

The data should be structured with one row per period per person. The data should be ordered so that for each person, the earliest period appears first in the data, the case period appears last.

There should be the same number of control periods and therefore rows of data for each person.

The following 3 variables are required:
Patient ID
Binary exposure indicator
Binary indicator for the outcome, 0 in control periods, 1 in case period for cases


More examples and an R package to come! Watch this space.