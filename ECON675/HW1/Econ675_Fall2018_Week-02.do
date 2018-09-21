********************************************************************************
** Econ-675: Applied Microeconometrics, Fall 2018
** Lecture 1: Causal Inference and Randomized Experiments
** Author: Matias D. Cattaneo
** Last update: 09-Sep-2018
********************************************************************************
clear all
set more off
set linesize 200

********************************************************************************
** Example: JTPA
********************************************************************************
use jtpa.dta, clear

********************************************************************************
** Convetional approach
********************************************************************************
** using TTEST
ttest earnings, by(treatment)
ttest earnings, by(treatment) unequal

** using stata
sum earnings if treatment==0
local N0 = r(N)
local mu0 = r(mean)
local sd0 = r(sd)
local V0 = r(Var)/r(N)

sum earnings if treatment==1
local N1 = r(N)
local mu1 = r(mean)
local sd1 = r(sd)
local V1 = r(Var)/r(N)

local tau = `mu1'-`mu0'
local v = sqrt(`V1'+`V0')
local T = `tau'/`v'
local pval = 2*normal(-abs(`T'))

local mu0 = round(`mu0', .01)
local mu1 = round(`mu1', .0001)
local sd0 = round(`sd0', .01)
local sd1 = round(`sd1', .0001)

di "Group Means (Std.Dev.)"
di "Control: `mu0' (`sd0') -- N=`N0'" 
di "Treatme: `mu1' (`sd1') -- N=`N1'" 
di "Test = `T' -- p-val = `pval'" 

* Compare to:
ttest earnings, by(treatment) unequal
ret list

** using mata
mata:
	y = st_data(., "earnings")
	t = st_data(., "treatment")
	N1 = sum(t:==1); N0 = rows(t)-N1
	mu1 = mean(select(y,t:==1)); mu1
	mu0 = mean(select(y,t:==0)); mu0
	v1 = variance(select(y,t:==1))/N1;
	v0 = variance(select(y,t:==0))/N0;
	
	T = (mu1-mu0)/sqrt(v1+v0)
	pval = 2*normal(-abs(T))
	
	st_numscalar("T", T); st_numscalar("pval", pval)
end

display "t-test = " T
display "pval   = " pval

********************************************************************************
** Randomization Inference approach
********************************************************************************
clear all
set obs 15
gen T = (runiform()>=.6)
mata: comb(15,sum(st_data(.,"T")))

gen Y = rnormal() + T*1
** using stata, ttest
ttest Y, by(T)
** using stata, permute
permute T diffmean=(r(mu_2)-r(mu_1)), reps(1000) nowarn: ttest Y, by(T) 
matrix pval = r(p)
display "p-val = " pval[1,1] 

** using mata, by hand
** careful: cvpermute() not random!
mata:
	y = st_data(., "Y"); t = st_data(., "T")
	
	totperm = comb(rows(t),sum(t)); totperm

	T = abs(mean(select(y,t:==1))-mean(select(y,t:==0)))

	S=3000; Tdist = J(min((S,totperm)),1,.)
	info = cvpermutesetup(t); s=1
	while ((t=cvpermute(info)) != J(0,1,.) & s<=S){
		Tdist[s,] = abs(mean(select(y,t:==1))-mean(select(y,t:==0)))
		s++
	}
	pval = sum(Tdist:>=T)/min((S,totperm)); pval
end

** Example: JTPA
use jtpa.dta, clear
permute treatment diffmean=(r(mu_2)-r(mu_1)), reps(500) nowarn: ttest earnings, by(treatment) 

* Compare to:
ttest earnings, by(treatment) unequal

********************************************************************************
** Experimental Design
********************************************************************************
local a = 0.05
local N1 = 100
local N2 = 500
local N3 = 1000
local se1 = sqrt(5 / `N1')
local se2 = sqrt(5 / `N2')
local se3 = sqrt(5 / `N3')

** Plot power functions
twoway (function y = 1 - normal(x/`se1'+invnormal(1-`a'/2)) + normal(x/`se1'-invnormal(1-`a'/2)), range(-1 1)) ///
	   (function y = 1 - normal(x/`se2'+invnormal(1-`a'/2)) + normal(x/`se2'-invnormal(1-`a'/2)), range(-1 1)) ///
	   (function y = 1 - normal(x/`se3'+invnormal(1-`a'/2)) + normal(x/`se3'-invnormal(1-`a'/2)), range(-1 1)), ///
	   legend(col(1) lab(1 "N=`N1'") lab(2 "N=`N2'") lab(3 "N=`N3'")) ///
	   yline(`a', lpattern(dash))

** Using ADO file
mypower earnings treatment

mypower earnings treatment, tau(500) rangeL(-1000) rangeU(1000)

















