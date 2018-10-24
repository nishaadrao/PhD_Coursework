## ECON675: ASSIGNMENT 3
## Q1: NLS
## Anirudh Yadav 
## 10/24/2018

######################################################################
# Load packages, clear workspace
######################################################################
rm(list = ls())             #clear workspace
library(foreach)            #for looping
library(data.table)         #for data manipulation
library(Matrix)             #fast matrix calcs
library(sandwich)           #for variance-covariance estimation 
library(xtable)             #for latex tables
options(scipen = 999)       #forces R to use normal numbers instead of scientific notation


######################################################################
# Input data, create additional covariates
######################################################################

# Get Piso Firme data
data <- as.data.table(read.csv('PhD_Coursework/ECON675/HW3/pisofirme.csv'))

# Create dependent variable for logistic regression
data[,s:= 1-dmissing]

# Create income regressor
data[,log_inc:= log(S_incomepc+1)]

######################################################################
# Q9(a): Estimate logistic regression
######################################################################

# Estimate model
mylogit <- glm(s ~ S_age + S_HHpeople + log_inc, data = data, family = "binomial")

b.hat   <- as.data.table(mylogit["coefficients"])

# Get robust standard errors
V.hat       <- vcovHC(mylogit, type = "HC1")
se.hat      <- as.data.table(sqrt(diag(V.hat)))

# Compute t-stats
t.stats       <- b.hat/se.hat

# Compute p-vals
n = nrow(data)
d = 4
p = round(2*pt(abs(t.stats[[1]]),df=n-d,lower.tail=FALSE),3)

# Compute CIs
CI.lower = b.hat - qnorm(0.975)*se.hat
CI.upper = b.hat + qnorm(0.975)*se.hat

# Mash results together
results  = as.data.frame(cbind(b.hat,se.hat,t.stats,p,CI.lower,CI.upper))
colnames(results) = c("Coef.","Std. Err.","t-stat","p-val","CI.lower","CI.upper")
rownames(results) = c("Const.", "S_age","S_HHpeople","log_inc")

# Get latex table output
xtable(results,digits=3)


