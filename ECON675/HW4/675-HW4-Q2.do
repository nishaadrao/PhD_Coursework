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
* [1] Difference in means
********************************************************************************

* Lalonde control
reg re78 treat if treat==1|treat==0 , hc2

* PSID control
reg re78 treat if treat==1|treat==2 , hc2

********************************************************************************
* [2] OLS
********************************************************************************

* Covariates A, Lalonde control
reg re78 treat age educ black hisp married nodegr log_re74 log_re75 if treat==1|treat==0 , hc2

* Covariates B, Lalonde control
reg re78 treat age educ black hisp married nodegr log_re74 log_re75 age_sq educ_sq u74 u75 if treat==1|treat==0 , hc2

* Covariates C, Lalonde control
reg re78 treat age educ black hisp married nodegr log_re74 log_re75 age_sq educ_sq u74 u75 age_cu black_u74 educ_log_re74 if treat==1|treat==0 , hc2

* Covariates A, PSID
reg re78 treat age educ black hisp married nodegr log_re74 log_re75 if treat==1|treat==2 , hc2

* Covariates B, PSID
reg re78 treat age educ black hisp married nodegr log_re74 log_re75 age_sq educ_sq u74 u75 if treat==1|treat==2 , hc2

* Covariates C, PSID
reg re78 treat age educ black hisp married nodegr log_re74 log_re75 age_sq educ_sq u74 u75 age_cu black_u74 educ_log_re74 if treat==1|treat==2 , hc2


********************************************************************************
* [3] Regression Imputation 
********************************************************************************

* Covariates A, Lalonde control
teffects ra (re78 age educ black hisp married nodegr log_re74 log_re75) (treat) if treat==1|treat==0 , ate 
teffects ra (re78 age educ black hisp married nodegr log_re74 log_re75) (treat) if treat==1|treat==0 , atet

* Covariates B, Lalonde control
teffects ra (re78 age educ black hisp married nodegr log_re74 log_re75 age_sq educ_sq u74 u75) (treat) if treat==1|treat==0 , ate 
teffects ra (re78 age educ black hisp married nodegr log_re74 log_re75 age_sq educ_sq u74 u75) (treat) if treat==1|treat==0 , atet 

* Covariates C, Lalonde control
teffects ra (re78 age educ black hisp married nodegr log_re74 log_re75 age_sq educ_sq u74 u75 age_cu black_u74 educ_log_re74) (treat) if treat==1|treat==0 , ate 
teffects ra (re78 age educ black hisp married nodegr log_re74 log_re75 age_sq educ_sq u74 u75 age_cu black_u74 educ_log_re74) (treat) if treat==1|treat==0 , atet 


* Covariates A, PSID control
teffects ra (re78 age educ black hisp married nodegr log_re74 log_re75) (treat) if treat==1|treat==2 , ate 
teffects ra (re78 age educ black hisp married nodegr log_re74 log_re75) (treat) if treat==1|treat==2 , atet

* Covariates B, PSID control
teffects ra (re78 age educ black hisp married nodegr log_re74 log_re75 age_sq educ_sq u74 u75) (treat) if treat==1|treat==2 , ate 
teffects ra (re78 age educ black hisp married nodegr log_re74 log_re75 age_sq educ_sq u74 u75) (treat) if treat==1|treat==2 , atet 

* Covariates C, PSID control
teffects ra (re78 age educ black hisp married nodegr log_re74 log_re75 age_sq educ_sq u74 u75 age_cu black_u74 educ_log_re74) (treat) if treat==1|treat==2 , ate 
teffects ra (re78 age educ black hisp married nodegr log_re74 log_re75 age_sq educ_sq u74 u75 age_cu black_u74 educ_log_re74) (treat) if treat==1|treat==2 , atet 


********************************************************************************
* [4] IPW
********************************************************************************

* Covariates A, Lalonde control
teffects ipw (re78) (treat age educ black hisp married nodegr log_re74 log_re75, logit) if treat==1|treat==0 , ate 
teffects ipw (re78) (treat age educ black hisp married nodegr log_re74 log_re75, logit) if treat==1|treat==0 , atet 

* Covariates B, Lalonde control
teffects ipw (re78) (treat age educ black hisp married nodegr log_re74 log_re75 age_sq educ_sq u74 u75, logit) if treat==1|treat==0 , ate 
teffects ipw (re78) (treat age educ black hisp married nodegr log_re74 log_re75 age_sq educ_sq u74 u75, logit) if treat==1|treat==0 , atet

* Covariates C, Lalonde control
teffects ipw (re78) (treat age educ black hisp married nodegr log_re74 log_re75 age_sq educ_sq u74 u75 age_cu black_u74 educ_log_re74, logit) if treat==1|treat==0 , ate 
teffects ipw (re78) (treat age educ black hisp married nodegr log_re74 log_re75 age_sq educ_sq u74 u75 age_cu black_u74 educ_log_re74, logit) if treat==1|treat==0 , atet 

* Covariates A, PSID control
teffects ipw (re78) (treat age educ black hisp married nodegr log_re74 log_re75, logit) if treat==1|treat==2 , ate 
teffects ipw (re78) (treat age educ black hisp married nodegr log_re74 log_re75, logit) if treat==1|treat==2 , atet 

* Covariates B, PSID control
teffects ipw (re78) (treat age educ black hisp married nodegr log_re74 log_re75 age_sq educ_sq u74 u75, logit) if treat==1|treat==2 , ate 
teffects ipw (re78) (treat age educ black hisp married nodegr log_re74 log_re75 age_sq educ_sq u74 u75, logit) if treat==1|treat==2 , atet

* Covariates C, PSID control
teffects ipw (re78) (treat age educ black hisp married nodegr log_re74 log_re75 age_sq educ_sq u74 u75 age_cu black_u74 educ_log_re74, logit) if treat==1|treat==2 , ate 
teffects ipw (re78) (treat age educ black hisp married nodegr log_re74 log_re75 age_sq educ_sq u74 u75 age_cu black_u74 educ_log_re74, logit) if treat==1|treat==2 , atet 

********************************************************************************
* [5] Doubly Robust
********************************************************************************

* Covariates A, Lalonde control
teffects ipwra (re78) (treat age educ black hisp married nodegr log_re74 log_re75, logit) if treat==1|treat==0 , ate 
teffects ipwra (re78) (treat age educ black hisp married nodegr log_re74 log_re75, logit) if treat==1|treat==0 , atet 

* Covariates B, Lalonde control
teffects ipwra (re78) (treat age educ black hisp married nodegr log_re74 log_re75 age_sq educ_sq u74 u75, logit) if treat==1|treat==0 , ate 
teffects ipwra (re78) (treat age educ black hisp married nodegr log_re74 log_re75 age_sq educ_sq u74 u75, logit) if treat==1|treat==0 , atet

* Covariates C, Lalonde control
teffects ipwra (re78) (treat age educ black hisp married nodegr log_re74 log_re75 age_sq educ_sq u74 u75 age_cu black_u74 educ_log_re74, logit) if treat==1|treat==0 , ate 
teffects ipwra (re78) (treat age educ black hisp married nodegr log_re74 log_re75 age_sq educ_sq u74 u75 age_cu black_u74 educ_log_re74, logit) if treat==1|treat==0 , atet 

* Covariates A, PSID control
teffects ipwra (re78) (treat age educ black hisp married nodegr log_re74 log_re75, logit) if treat==1|treat==2 , ate 
teffects ipwra (re78) (treat age educ black hisp married nodegr log_re74 log_re75, logit) if treat==1|treat==2 , atet 

* Covariates B, PSID control
teffects ipwra (re78) (treat age educ black hisp married nodegr log_re74 log_re75 age_sq educ_sq u74 u75, logit) if treat==1|treat==2 , ate 
teffects ipwra (re78) (treat age educ black hisp married nodegr log_re74 log_re75 age_sq educ_sq u74 u75, logit) if treat==1|treat==2 , atet

* Covariates C, PSID control
teffects ipwra (re78) (treat age educ black hisp married nodegr log_re74 log_re75 age_sq educ_sq u74 u75 age_cu black_u74 educ_log_re74, logit) if treat==1|treat==2 , ate 
teffects ipwra (re78) (treat age educ black hisp married nodegr log_re74 log_re75 age_sq educ_sq u74 u75 age_cu black_u74 educ_log_re74, logit) if treat==1|treat==2 , atet 

