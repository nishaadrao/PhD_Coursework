## ECON675: ASSIGNMENT 5
## Q2: WEAK INSTRUMENTS SIMULATIONS
## Anirudh Yadav 
## 11/19/2018

######################################################################
# Load packages, clear workspace
######################################################################
rm(list = ls())             #clear workspace
library(foreach)            #for looping
library(data.table)         #for data manipulation
library(Matrix)             #fast matrix calcs
library(ggplot2)            #for pretty plots
library(sandwich)           #for variance-covariance estimation 
library(xtable)             #for latex tables
library(boot)               #for bootstrapping
library(mvtnorm)            #for MVN stuff
library(AER)                #for IV regressions
options(scipen = 999)       #forces R to use normal numbers instead of scientific notation

######################################################################
# Generate random data for each simulation
######################################################################
N     = 200
M     = 5000
SIGMA = matrix(c(1,0,0,0,1,0.99,0,0.99,1),3,3)

set.seed(1234)

# Generate Z, U, V
W     = replicate(M,rmvnorm(N, mean = c(0,0,0), sigma = SIGMA, method="chol"))

# Get Y (assuming that beta=0, Y=U)
Y     = W[,2,]

# Generate X matrix for each value of gamma
gamma.vec = sqrt((1/N)*c(0,0.25,9,99))
X         = lapply(1:length(gamma.vec),function(i) gamma.vec[i]*W[,1,]+W[,3,])

######################################################################
# Compute OLS statistics for each gamma, and simulation
######################################################################

# Run OLS for each gamma and each simulation -- this spits out 5000 lm's for each gamma
ols.big = foreach(j=1:length(gamma.vec)) %do%
  lapply(1:M, function(i) lm(Y[,i]~X[[j]][,i]-1))

# Extract point estimates, standard errors, t-stats
ols.beta = foreach(j=1:length(gamma.vec)) %do%
  sapply(1:M, function(i) ols.big[[j]][[i]]$coefficients)

ols.se   = foreach(j=1:length(gamma.vec)) %do%
  sapply(1:M, function(i) coef(summary(ols.big[[j]][[i]]))[,"Std. Error"])

ols.t    = sapply(1:length(gamma.vec),function(j) ols.beta[[j]]/ols.se[[j]])

ols.rej  = ifelse(ols.t>1.96,1,0)

# Compute desired summary statistics across the simulations (spits out a list containing results for each gamma)
ols.results = foreach(j=1:length(gamma.vec)) %do%
  rbind(c(mean(ols.beta[[j]]),sd(ols.beta[[j]]),quantile(ols.beta[[j]], probs = c(0.1, 0.5 ,0.9))),
  c(mean(ols.se[[j]]),sd(ols.se[[j]]),quantile(ols.se[[j]], probs = c(0.1, 0.5 ,0.9))),
  c(mean(ols.rej[,j]),sd(ols.rej[,j]),quantile(ols.rej[,j], probs = c(0.1, 0.5 ,0.9))))

# Remove big objects!
rm(ols.big,ols.beta,ols.se,ols.t,ols.rej)

######################################################################
# Compute 2SLS statistics for each gamma, and simulation
######################################################################

# Run 2SLS for each gamma and each simulation -- this spits out 5000 ivreg's for each gamma
# WATCH OUT: this takes a minute or so!
iv.big = foreach(j=1:length(gamma.vec)) %do%
  lapply(1:M, function(i) ivreg(Y[,i]~X[[j]][,i]-1|W[,1,i]))

# Extract point estimates, standard errors, t-stats
iv.beta = foreach(j=1:length(gamma.vec)) %do%
  sapply(1:M, function(i) iv.big[[j]][[i]]$coefficients)

iv.se   = foreach(j=1:length(gamma.vec)) %do%
  sapply(1:M, function(i) summary(iv.big[[j]][[i]])[["coefficients"]][,"Std. Error"])

iv.t    = sapply(1:length(gamma.vec),function(j) iv.beta[[j]]/iv.se[[j]])

iv.rej  = ifelse(iv.t>1.96,1,0)

# Run first-stage regression and extract F-statistics
iv.f       = foreach(j=1:length(gamma.vec)) %do%
  sapply(1:M, function(i) summary(lm(X[[j]][,i]~W[,1,i]-1))$fstatistic[1])

# Combine results for each gamma
iv.results = foreach(j=1:length(gamma.vec)) %do%
  rbind(c(mean(iv.beta[[j]]),sd(iv.beta[[j]]),quantile(iv.beta[[j]], probs = c(0.1, 0.5 ,0.9))),
        c(mean(iv.se[[j]]),sd(iv.se[[j]]),quantile(iv.se[[j]], probs = c(0.1, 0.5 ,0.9))),
        c(mean(iv.rej[,j]),sd(iv.rej[,j]),quantile(iv.rej[,j], probs = c(0.1, 0.5 ,0.9))),
        c(mean(iv.f[[j]]),sd(iv.f[[j]]),quantile(iv.f[[j]], probs = c(0.1, 0.5 ,0.9))))

# Remove big objects!
rm(iv.big,iv.beta,iv.se,iv.t,iv.rej,iv.f)


