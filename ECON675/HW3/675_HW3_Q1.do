********************************************************************************
* ECON675: ASSIGNMENT 3
* Q1: NLS
* Anirudh Yadav
* 10/24/2018
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
import delimited using "$dir/PhD_Coursework/ECON675/HW3/pisofirme.csv"

* Generate additional variables
gen       s = 1-dmissing
gen log_inc = log(s_incomepc+1)


********************************************************************************
* Q9(a): Estimate logistic regression
********************************************************************************
logit s s_age s_hhpeople log_inc, robust
