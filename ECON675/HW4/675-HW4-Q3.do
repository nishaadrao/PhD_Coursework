********************************************************************************
* ECON675: ASSIGNMENT 4
* Q3: POST-MODEL SELECTION INFERENCE
* Anirudh Yadav
* 11/09/2018
********************************************************************************


********************************************************************************
* Preliminaries
********************************************************************************
clear all
set more off

* Set working directory 
global dir "/Users/Anirudh/Desktop/GitHub"


set seed 22
set obs 50

********************************************************************************
* [1] Summary stats and density plots
********************************************************************************

* number of replications
local M = 1000
set matsize 11000

* empty matrices to store estimates and indicator of coverage
matrix est = J(`M',3,.)
matrix cov = J(`M',3,.)

* initial values we will replace during replication
gen x = rnormal(0,1) 
gen z = .85*x + sqrt(1-.85)*rnormal(0,1)
gen eps = rnormal(0,1)
gen y = 1 + .5*x + z + eps

* loop for M replications
forvalues i = 1/`M'{
	qui replace x = rnormal(0,1) 
	qui replace z = .85*x + sqrt(1-.85)*rnormal(0,1)
	qui replace eps = rnormal(0,1)
	qui replace y = 1 + .5*x + z + eps
	
	* long regression
	qui reg y x z, r
	
	* extract first estimate
	local beta_hat = _b["x"]
	matrix est[`i',1] = `beta_hat'
	
	* get SE and calculate coverage of true beta_0 = .5
	local se_hat = _se["x"]
	local lb_hat = `beta_hat' - 1.96 * `se_hat'
	local ub_hat = `beta_hat' + 1.96 * `se_hat'
	local cov_hat = (.5 >= `lb_hat') & (.5 <= `ub_hat')
	matrix cov[`i',1] = `cov_hat'
	
	* save gamma over se gamma
	local gamma_hat = _b["z"]
	local gamma_se  = _se["z"]
	local tstat = `gamma_hat'/`gamma_se'
		
	* short regression
	qui reg y x, r
	local beta_tilde = _b["x"]
	matrix est[`i',2] = `beta_tilde'
	
	* get SE and calculate coverage of true beta_0 = .5
	local se_tilde = _se["x"]
	local lb_tilde = `beta_tilde' - 1.96 * `se_tilde'
	local ub_tilde = `beta_tilde' + 1.96 * `se_tilde'
	local cov_tilde = (.5 >= `lb_tilde') & (.5 <= `ub_tilde')
	matrix cov[`i',2] = `cov_tilde'
	
	* third estimate
	local beta_check = cond(`tstat' >= 1.96, `beta_hat', `beta_tilde')	
	matrix est[`i',3] = cond(`tstat' >= 1.96, `beta_hat', `beta_tilde')	
	matrix cov[`i',3] = cond(`tstat' >= 1.96, `cov_hat', `cov_tilde') 
}

* turn results into variables
svmat est
svmat cov

* drop old data
drop x
drop z
drop eps
drop y

* rename variables
rename est1 beta_hat
rename est2 beta_tilde
rename est3 beta_check
rename cov1 cov_hat
rename cov2 cov_tilde
rename cov3 cov_check

* write summary statistics to latex
outreg2 using q3.tex, replace sum(log) ///
	keep(beta_hat beta_tilde beta_check) ///
	eqkeep(min mean median max) ///
	dec(2)

* kernel densities
twoway kdensity beta_hat, k(epanechnikov) || ///
 kdensity beta_tilde, k(epanechnikov) || ///
 kdensity beta_check, k(epanechnikov) ///
 leg(lab(1 "beta_hat") lab(2 "beta_tilde") lab(3 "beta_check")) ///
 ytitle("Density") xtitle("")
	 

********************************************************************************
* [2] Coverage rates
********************************************************************************

* calculate these here, report them in LaTeX
sum(cov_hat)
sum(cov_tilde)
sum(cov_check)
