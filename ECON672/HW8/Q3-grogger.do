********************************************************************************
** ECON672 PS8 Question 3 - Probability of Arrest
*
* Anirudh Yadav - 23 March 2018
*
*
********************************************************************************

clear all
use "/Users/Anirudh/Desktop/PhD/ECON672/statafiles-wooldridge/grogger.dta"
*ssc install estout, replace

********************************************************************************
* Part (a)
********************************************************************************

gen arr86 = 0
replace arr86 = 1 if narr86>0

eststo: quietly reg arr86 pcnv avgsen tottime ptime86 inc86 black hispan born60, r

********************************************************************************
* Part (b) 
********************************************************************************

* Non-robust test
reg arr86 pcnv avgsen tottime ptime86 inc86 black hispan born60
test avgsen tottime
* Robust test
reg arr86 pcnv avgsen tottime ptime86 inc86 black hispan born60, r
test avgsen tottime

********************************************************************************
* Part (c)
********************************************************************************
eststo: quietly probit arr86 c.pcnv c.avgsen c.tottime c.ptime86 c.inc86 i.black i.hispan i.born60, r

* Compute estimated effect at given margins of increasing pcnv from 0.25 to 0.75
margins, at((mean) _all black=1 hispan=0 born60=1 pcnv=(0.25 0.75)) 

* Compute difference in estimate probabilities
matrix b = r(b)
display b[1,2] - b[1,1]



********************************************************************************
* Part (d)
********************************************************************************

eststo: quietly probit arr86 pcnv avgsen tottime ptime86 inc86 black hispan born60 pcnvsq pt86sq inc86sq, r
test pcnvsq pt86sq inc86sq


** LR TEST ** NOTE: the LR test doesn't work with robust SEs...
* Unconstrained model
quietly probit arr86 pcnv avgsen tottime ptime86 inc86 black hispan born60 pcnvsq pt86sq inc86sq
estimates store full

* Constrained model
quietly probit arr86 pcnv avgsen tottime ptime86 inc86 black hispan born60

* Test
lrtest full


* Export regression results
esttab using "/Users/Anirudh/Desktop/PhD/ECON672/Problem Sets/PS8/TeX/Results.tex", se replace

