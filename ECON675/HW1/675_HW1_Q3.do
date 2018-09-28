********************************************************************************
* ECON675: ASSIGNMENT 1
* Q3: ANALYSIS OF EXPERIMENTS
* Anirudh Yadav
* 8/26/2018
********************************************************************************

********************************************************************************
* Preliminaries
********************************************************************************
clear all
set more off

* Set working directory 
global dir "/Users/Anirudh/Desktop/GitHub"


********************************************************************************
* Import data
********************************************************************************

* Import LaLonde data
import delimited using "$dir/PhD_Coursework/ECON675/HW1/LaLonde_1986.csv"


********************************************************************************
* Q1(a): Difference in means estimator
********************************************************************************

sum earn78 if treat==0
local N0 = r(N)
local mu0 = r(mean)
local sd0 = r(sd)
local V0 = r(Var)/r(N)

sum earn78 if treat==1
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

di "`tau'"

********************************************************************************
* Q1(b): Conservative confidence intervals
********************************************************************************

local CIlower = `tau' - invnormal(0.975)*`v'
local CIupper = `tau' + invnormal(0.975)*`v'

di "`CIlower'"
di "`CIupper'"


********************************************************************************
* Q2(a): Fisher exact p-values
********************************************************************************

* Using difference in means estimator
permute treat diffmean=(r(mu_2)-r(mu_1)), reps(1999) nowarn: ttest earn78, by(treat) 
matrix pval = r(p)
display "p-val = " pval[1,1]

* Using KS statistic
permute treat ks=r(D), reps(1999) nowarn: ksmirnov earn78, by(treat) 
matrix pval = r(p)
display "p-val = " pval[1,1]


********************************************************************************
* Q2(b): Fisher confidence intervals
********************************************************************************

* Infer missing values under the null of constant treatment effect
gen     Y1_imputed = earn78
replace Y1_imputed = earn78 + `tau' if treat==0

gen     Y0_imputed = earn78
replace Y0_imputed = earn78 - `tau' if treat==1

* Write program to put into bootstrap function
program define meandiff, rclass
	summarize   Y1_imputed if treat==1
	local 		tau1 = r(mean)
	sum 		Y0_imputed if treat==0
	local 		tau0 = r(mean)
	return      scalar meandiff = `tau1' - `tau0'
end

* Run bootstrap function using meandiff program
bootstrap diff = r(meandiff), reps(1999): meandiff

********************************************************************************
* Q3(a): Plot power function
********************************************************************************
local alpha = 0.05
local z     = invnormal(1-`alpha'/2)

* Plot power function
twoway (function y = 1 - normal(`z'- x/`v') + normal(-`z'-x/`v'), range(-4000 //
		4000)), ///
	   yline(`alpha', lpattern(dash)) yti("Power") xti("tau")
	   
********************************************************************************
* Q3(b): Sample size calculation 
********************************************************************************
mata:

mata clear
 
  function mypowerfunc(N, S0, S1, p, tau){

  return(1 - normal(invnormal(0.975)-tau/sqrt(1/N*S1*(1/p)+1/N*S0*(1/(1-p)))) +
  normal(-invnormal(0.975)-tau/sqrt(1/N*S1*(1/p)+1/N*S0*(1/(1-p)))) -0.8)
 
		}
		
	y = st_data(., "earn78")
	t = st_data(., "treat")	
		
	S1 = variance(select(y,t:==1))
	S0 = variance(select(y,t:==0))
    p     = 2/3
	tau   = 1000

 
	mm_root(x=., &mypowerfunc(), 1000, 1500, 0, 10000, S0,S1, p ,tau)
      
	x
	  
end 
 
