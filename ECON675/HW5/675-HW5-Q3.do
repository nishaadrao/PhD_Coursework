********************************************************************************
* ECON675: ASSIGNMENT 5
* Q3: WEAK INSTRUMENTS -- EMPIRICAL STUDIES
* Anirudh Yadav
* 11/19/2018
********************************************************************************


********************************************************************************
* Preliminaries
********************************************************************************
clear all
set more off

* Set working directory 
global dir "/Users/Anirudh/Desktop/GitHub"


********************************************************************************
* Import AK data
********************************************************************************

use "$dir/PhD_Coursework/ECON675/HW5/Angrist_Krueger.dta"

********************************************************************************
* [3.1] Run AK regressions
********************************************************************************
eststo ols1: reg l_w_wage educ non_white married SMSA i.region i.YoB_ld, r
eststo ols2: reg l_w_wage educ non_white married SMSA i.region i.YoB_ld age_q age_sq, r
eststo  iv1: ivregress 2sls l_w_wage non_white married SMSA i.region i.YoB_ld (educ = i.QoB##i.YoB_ld),r
eststo  iv2: ivregress 2sls l_w_wage non_white married SMSA i.region i.YoB_ld age_q age_sq (educ = i.QoB##i.YoB_ld), r

esttab ols1 ols2 iv1 iv2 using "$dir/PhD_Coursework/ECON675/HW5/q3_ak_results.tex", keep(educ non_white SMSA married age_q age_sq) se nostar

********************************************************************************
* [3.2] Run BJB permutation regressions
********************************************************************************
capture program drop IV_quick
program define IV_quick, rclass
    syntax varlist(max=1) [, model(integer 1) ]
	local x "`varlist'"
	
	if (`model' == 1) {
	    capture drop educ_hat
		qui reg educ non_white married SMSA i.region i.YoB_ld i.YoB_ld##i.`x'
		predict educ_hat
		qui reg l_w_wage educ_hat non_white married SMSA i.region i.YoB_ld
		return scalar beta = _b[educ_hat]
	}
	if (`model' == 2) {
	    capture drop educ_hat
		qui reg educ non_white married SMSA age_q age_sq i.region i.YoB_ld i.YoB_ld##i.`x'
		predict educ_hat
		qui reg l_w_wage educ_hat non_white married SMSA age_q age_sq i.region i.YoB_ld
		return scalar beta = _b[educ_hat]
	}
end


permute QoB TSLS_1_b = r(beta), reps(500) seed(123) saving("$dir/PhD_Coursework/ECON675/HW5/premute1.dta", replace): ///
    IV_quick QoB, model(1)
	
permute QoB TSLS_2_b = _b[educ], reps(500) seed(123) saving("$dir/PhD_Coursework/ECON675/HW5/premute2.dta", replace): ///
    ivregress 2sls l_w_wage non_white married SMSA age_q age_sq i.region i.YoB_ld ///
    (educ = i.YoB_ld##i.QoB)

clear all
use "$dir/PhD_Coursework/ECON675/HW5/premute1.dta"
sum TSLS_1_b

clear all
use "$dir/PhD_Coursework/ECON675/HW5/premute1.dta"
sum TSLS_2_b
