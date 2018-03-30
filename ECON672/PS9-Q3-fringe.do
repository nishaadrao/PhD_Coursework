********************************************************************************
** ECON672 PS9 Question 3 - Fringe Benefits
*
* Anirudh Yadav - 30 March 2018
*
*
********************************************************************************

clear all
use "/Users/Anirudh/Desktop/PhD/ECON672/statafiles-wooldridge/fringe.dta"
*ssc install estout, replace

********************************************************************************
* Part (a)
********************************************************************************

eststo: quietly reg hrbens exper age educ tenure married male white nrtheast nrthcen south union


********************************************************************************
* Part (b)
********************************************************************************

eststo: quietly tobit hrbens exper age educ tenure married male white nrtheast nrthcen south union, ll

********************************************************************************
* Part (c)
********************************************************************************

eststo: quietly tobit hrbens exper age educ tenure married male white nrtheast nrthcen south union expersq tenuresq, ll

** Wald (?) test **
test expersq tenuresq

** LR TEST ** this should be the preferred test I think, but can't use robust SEs
* Unconstrained model
quietly tobit hrbens exper age educ tenure married male white nrtheast nrthcen south union expersq tenuresq, ll
estimates store full

* Constrained model
quietly tobit hrbens exper age educ tenure married male white nrtheast nrthcen south union, ll

* Test
lrtest full

esttab using "/Users/Anirudh/Desktop/PhD/ECON672/Problem Sets/PS9/TeX/Results.tex", se replace



