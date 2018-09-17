//ECON675: ASSIGNMENT 1
//Q2: IMPLEMENTING LEAST-SQUARES ESTIMATORS
//Anirudh Yadav
// 8/16/2018

//------ Data import and processing ------//
import delimited /Users/Anirudh/Desktop/GitHub/PhD_Coursework/ECON675/HW1/LaLonde_1986.csv,clear

gen educsq=educ^2
gen black_earn74 = black*earn74


//------ Run OLS regression -------------//
reg earn78 treat black age educ educsq earn74 black_earn74 u74 u75, r

