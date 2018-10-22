********************************************************************************
* ECON641: PS1 
* Q8: GRAVITY
* Anirudh Yadav
* 10/22/2018
********************************************************************************


********************************************************************************
* Preliminaries
********************************************************************************
clear all
set more off

* Set working directory 
global dir "/Users/Anirudh/Desktop/GitHub//PhD_Coursework/ECON641/HW1"

********************************************************************************
* Import data create additional covariates
********************************************************************************

* Import gravity dataset
use "$dir/col_regfile09.dta"

* Create depvar as in AvW
gen   y   = flow/(gdp_o*gdp_d)
gen log_y = ln(y)

* Log distance
gen log_d = ln(distw)

* No common language
gen nocomml = 1-comlang_off

* No colonial ties
gen nocol = 1-col_hist

* No contiguity
gen nocontig = 1-contig

* Importer fixed effect
encode iso_d, gen(m_i)

* Exporter fixed effect
encode iso_o, gen(m_e)


********************************************************************************
* Run OLS for year 2000
********************************************************************************
eststo: reg log_y log_d nocontig nocomml nocol i.m_i i.m_e if log_y!=. & year==2000, robust


********************************************************************************
* Run PPML for year 2000
********************************************************************************

* Positive trade sample
eststo: poisson y log_d nocontig nocomml nocol i.m_i i.m_e if flow>0 & year==2000, robust


* Incl. zeros
eststo: poisson y log_d nocontig nocomml nocol i.m_i i.m_e if flow!=. & year==2000, robust


********************************************************************************
* Make table for tex file
********************************************************************************
esttab using Q8_results.tex, se nostar drop(*m_i *m_e) label replace booktabs
