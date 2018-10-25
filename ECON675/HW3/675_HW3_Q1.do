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


********************************************************************************
* Q9(b): Bootstrap statistics
********************************************************************************
logit s s_age s_hhpeople log_inc, vce(bootstrap, reps(999))


********************************************************************************
* Q9(b): Predicted probabilities
********************************************************************************
logit s s_age s_hhpeople log_inc, vce(robust)

* predict propensity score
predict p

* plot kernel density estimates
twoway histogram p || kdensity p, k(gaussian) || ///
 kdensity p, k(epanechnikov) || kdensity p, k(triangle) ///
 leg(lab(1 "Propensity Score") lab(2 "Gaussian") ///
	 lab(3 "Epanechnikov") lab(4 "Triangle"))
