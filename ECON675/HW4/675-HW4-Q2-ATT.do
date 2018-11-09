********************************************************************************
* ECON675: ASSIGNMENT 4
* Q2: ESTIMATING AVERAGE TREATMENT EFFECTS
* Anirudh Yadav
* 11/07/2018
********************************************************************************


********************************************************************************
* Preliminaries
********************************************************************************
clear all
set more off

* Set working directory 
global dir "/Users/Anirudh/Desktop/GitHub"


********************************************************************************
* Import data, create additional covariates
********************************************************************************

* Import LaLonde data
import delimited using "$dir/PhD_Coursework/ECON675/HW4/LaLonde_all.csv"

* Generate additional covariates 
gen log_re74 = log(re74+1)
gen log_re75 = log(re75+1)
gen age_sq   = age^2
gen age_cu   = age^3
gen educ_sq  = educ^2
gen black_u74 = black*u74
gen educ_log_re74 = educ*log_re74


********************************************************************************
* [3] Regression Imputation 
********************************************************************************

* Covariates A, Lalonde control
eststo r1: teffects ra (re78 age educ black hisp married nodegr log_re74 log_re75) (treat) if treat==1|treat==0 , atet

* Covariates B, Lalonde control
eststo r2: teffects ra (re78 age educ black hisp married nodegr log_re74 log_re75 age_sq educ_sq u74 u75) (treat) if treat==1|treat==0 , atet 

* Covariates C, Lalonde control
eststo r3: teffects ra (re78 age educ black hisp married nodegr log_re74 log_re75 age_sq educ_sq u74 u75 age_cu black_u74 educ_log_re74) (treat) if treat==1|treat==0 , atet 

esttab r1 r2 r3 using Q2_att.csv, se nostar keep(r1vs0.treat) wide noparentheses nonumber noobs plain nomtitles replace


* Covariates A, PSID control
eststo r4: teffects ra (re78 age educ black hisp married nodegr log_re74 log_re75) (treat) if treat==1|treat==2 , atet

* Covariates B, PSID control
eststo r5: teffects ra (re78 age educ black hisp married nodegr log_re74 log_re75 age_sq educ_sq u74 u75) (treat) if treat==1|treat==2 , atet 

* Covariates C, PSID control
eststo r6: teffects ra (re78 age educ black hisp married nodegr log_re74 log_re75 age_sq educ_sq u74 u75 age_cu black_u74 educ_log_re74) (treat) if treat==1|treat==2 , atet 


esttab r4 r5 r6 using Q2_att.csv, se nostar keep(r2vs1.treat) wide noparentheses nonumber noobs plain nomtitles append

********************************************************************************
* [4] IPW
********************************************************************************

* Covariates A, Lalonde control
eststo i1: teffects ipw (re78) (treat age educ black hisp married nodegr log_re74 log_re75, logit) if treat==1|treat==0 , atet 

* Covariates B, Lalonde control
eststo i2: teffects ipw (re78) (treat age educ black hisp married nodegr log_re74 log_re75 age_sq educ_sq u74 u75, logit) if treat==1|treat==0 , atet

* Covariates C, Lalonde control
eststo i3: teffects ipw (re78) (treat age educ black hisp married nodegr log_re74 log_re75 age_sq educ_sq u74 u75 age_cu black_u74 educ_log_re74, logit) if treat==1|treat==0 , atet 

esttab i1 i2 i3 using Q2_att.csv, se nostar keep(r1vs0.treat) wide noparentheses nonumber noobs plain nomtitles append


* Covariates A, PSID control [doesn't converge, so set maxiter = 50!!!]
eststo i4: teffects ipw (re78) (treat age educ black hisp married nodegr log_re74 log_re75, logit) if treat==1|treat==2 , atet iterate(50)

* Covariates B, PSID control [first need to drop obs with very low prop scores]
teffects ipw (re78) (treat age educ black hisp married nodegr log_re74 log_re75 age_sq educ_sq u74 u75, logit) if treat==1|treat==2 , ate osample(viol)
eststo i5: teffects ipw (re78) (treat age educ black hisp married nodegr log_re74 log_re75 age_sq educ_sq u74 u75, logit) if treat==1|treat==2 & viol==0 , atet iter(50)

* Covariates C, PSID control
eststo i6: teffects ipw (re78) (treat age educ black hisp married nodegr log_re74 log_re75 age_sq educ_sq u74 u75 age_cu black_u74 educ_log_re74, logit) if treat==1|treat==2 , atet osample(violl)
eststo i6: teffects ipw (re78) (treat age educ black hisp married nodegr log_re74 log_re75 age_sq educ_sq u74 u75 age_cu black_u74 educ_log_re74, logit) if treat==1|treat==2 & violl==0, atet iter(50) 


esttab i4 i5 i6 using Q2_att.csv, se nostar keep(r2vs1.treat) wide noparentheses nonumber noobs plain nomtitles append

********************************************************************************
* [5] Doubly Robust
********************************************************************************

* Covariates A, Lalonde control
eststo dr1: teffects ipwra (re78) (treat age educ black hisp married nodegr log_re74 log_re75, logit) if treat==1|treat==0 , atet 

* Covariates B, Lalonde control
eststo dr2: teffects ipwra (re78) (treat age educ black hisp married nodegr log_re74 log_re75 age_sq educ_sq u74 u75, logit) if treat==1|treat==0 , atet

* Covariates C, Lalonde control
eststo dr3: teffects ipwra (re78) (treat age educ black hisp married nodegr log_re74 log_re75 age_sq educ_sq u74 u75 age_cu black_u74 educ_log_re74, logit) if treat==1|treat==0 , atet 

esttab dr1 dr2 dr3 using Q2_att.csv, se nostar keep(r1vs0.treat) wide noparentheses nonumber noobs plain nomtitles append


* Covariates A, PSID control
eststo dr4: teffects ipwra (re78) (treat age educ black hisp married nodegr log_re74 log_re75, logit) if treat==1|treat==2 , atet iter(25) 

* Covariates B, PSID control
eststo dr5: teffects ipwra (re78) (treat age educ black hisp married nodegr log_re74 log_re75 age_sq educ_sq u74 u75, logit) if treat==1|treat==2 & viol==0 , atet iter(25)

* Covariates C, PSID control
eststo dr6: teffects ipwra (re78) (treat age educ black hisp married nodegr log_re74 log_re75 age_sq educ_sq u74 u75 age_cu black_u74 educ_log_re74, logit) if treat==1|treat==2 & violl==0 , atet iter(25)

esttab dr4 dr5 dr6 using Q2_att.csv, se nostar keep(r2vs1.treat) wide noparentheses nonumber noobs plain nomtitles append

********************************************************************************
* [6] Nearest Neighbour Matching
********************************************************************************

* Covariates A, Lalonde control
eststo nn1: teffects nnmatch (re78 age educ black hisp married nodegr log_re74 log_re75) (treat) if treat==1|treat==0 , atet nneighbor(1) metric(maha)

* Covariates B, Lalonde control
eststo nn2: teffects nnmatch (re78 age educ black hisp married nodegr log_re74 log_re75 age_sq educ_sq u74 u75) (treat) if treat==1|treat==0 , atet nneighbor(1) metric(maha)

* Covariates C, Lalonde control
eststo nn3: teffects nnmatch (re78 age educ black hisp married nodegr log_re74 log_re75 age_sq educ_sq u74 u75 age_cu black_u74 educ_log_re74) (treat) if treat==1|treat==0 , atet nneighbor(1) metric(maha)

esttab nn1 nn2 nn3 using Q2_att.csv, se nostar keep(r1vs0.treat) wide noparentheses nonumber noobs plain nomtitles append


* Covariates A, PSID control
eststo nn4: teffects nnmatch (re78 age educ black hisp married nodegr log_re74 log_re75) (treat) if treat==1|treat==2 , atet nneighbor(1) metric(maha)

* Covariates B, PSID control
eststo nn5: teffects nnmatch (re78 age educ black hisp married nodegr log_re74 log_re75 age_sq educ_sq u74 u75) (treat) if treat==1|treat==2 , atet nneighbor(1) metric(maha)

* Covariates C, PSID control
eststo nn6: teffects nnmatch (re78 age educ black hisp married nodegr log_re74 log_re75 age_sq educ_sq u74 u75 age_cu black_u74 educ_log_re74) (treat) if treat==1|treat==2 , atet nneighbor(1) metric(maha)


esttab nn4 nn5 nn6 using Q2_att.csv, se nostar keep(r2vs1.treat) wide noparentheses nonumber noobs plain nomtitles append

********************************************************************************
* [7] PS matching
********************************************************************************

* Covariates A, Lalonde control
eststo ps1: teffects psmatch (re78) (treat age educ black hisp married nodegr log_re74 log_re75, logit) if treat==1|treat==0 , atet 

* Covariates B, Lalonde control
eststo ps2: teffects psmatch (re78) (treat age educ black hisp married nodegr log_re74 log_re75 age_sq educ_sq u74 u75, logit) if treat==1|treat==0 , atet

* Covariates C, Lalonde control
eststo ps3: teffects psmatch (re78) (treat age educ black hisp married nodegr log_re74 log_re75 age_sq educ_sq u74 u75 age_cu black_u74 educ_log_re74, logit) if treat==1|treat==0 , atet 

esttab ps1 ps2 ps3 using Q2_att.csv, se nostar keep(r1vs0.treat) wide noparentheses nonumber noobs plain nomtitles append


* Covariates A, PSID control
eststo ps4: teffects psmatch (re78) (treat age educ black hisp married nodegr log_re74 log_re75, logit) if treat==1|treat==2 , atet 

* For the PSID samples below there are some prop scores too close to 1.
* First I need to run the treatment models, identify the respondents w/ problematic prop scores -- this will cause the code to break
* Then I drop the violators and estimate the treatment effects
teffects psmatch (re78) (treat age educ black hisp married nodegr log_re74 log_re75 age_sq educ_sq u74 u75 age_cu black_u74 educ_log_re74, logit) if treat==1|treat==2 , ate osample(viol3)


* Covariates B, PSID control
eststo ps5: teffects psmatch (re78) (treat age educ black hisp married nodegr log_re74 log_re75 age_sq educ_sq u74 u75, logit) if treat==1|treat==2 & viol==0, atet 

* Covariates C, PSID control
eststo ps6: teffects psmatch (re78) (treat age educ black hisp married nodegr log_re74 log_re75 age_sq educ_sq u74 u75 age_cu black_u74 educ_log_re74, logit) if treat==1|treat==2 & violl==0 , atet 

esttab ps4 ps5 ps6 using Q2_att.csv, se nostar keep(r2vs1.treat) wide noparentheses nonumber noobs plain nomtitles append
