********************************************************************************
* ECON675: ASSIGNMENT 1
* Q2: IMPLEMENTING LEAST-SQUARES ESTIMATORS
* Anirudh Yadav
* 8/16/2018
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
import delimited using "$dir/PhD_Coursework/ECON675/HW1/LaLonde_1986.csv"

* Generate additional covariates
gen educsq=educ^2
gen black_earn74 = black*earn74
gen ones = 1

********************************************************************************
* Q4: Matrix implementation of OLS
********************************************************************************
mata:

// Create data matricies
X 		= st_data(.,("ones", "treat", "black", "age", "educ", "educsq", "earn74","black_earn74", "u74","u75"))
Y 		= st_data(.,("earn78"))

// Compute OLS point estimator
M 		= invsym(cross(X,X))
betahat = M*cross(X,Y)

// Construct diagonal matrix of squared residuals
U 		= Y - X*betahat
D		= diag(U*U')

// Compute asymptotic White var-cov matrix 
n 		= rows(X)
d 		= cols(X)

V 		= n*M*X'*D*X*M

// Compute standard errors
se 		= sqrt(diag(V)/n)
se  	= diagonal(se)

// Compute t-statistics (element-wise division)
t       = betahat:/se

// Compute p-values
p       = 2*ttail(n-d, t)

// Compute 95% confidence intervals
CIlower = betahat - invnormal(0.975)*se
CIupper = betahat + invnormal(0.975)*se

betahat, se, t , p , CIlower, CIupper
end

********************************************************************************
* Q5(b): compute OLS results using reg function
********************************************************************************

reg earn78 treat black age educ educsq earn74 black_earn74 u74 u75, r

* NOTE that the differences in se's is because the "r" option implements the
* d.f adjustment; i.e. se(reg) = n/(n-d)*se(mata). 
* I could easily implement the d.f. adjustment in the mata implementation, but
* I think it's nice to see the comparison.

