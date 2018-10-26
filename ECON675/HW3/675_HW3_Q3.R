## ECON675: ASSIGNMENT 3
## Q3: WHEN BOOTSTRAP FAILS
## Anirudh Yadav 
## 10/26/2018

######################################################################
# Load packages, clear workspace
######################################################################
rm(list = ls())             #clear workspace
library(foreach)            #for looping
library(data.table)         #for data manipulation
library(dplyr)              #melting for ggplot
library(Matrix)             #fast matrix calcs
library(ggplot2)            #for pretty plots
library(sandwich)           #for variance-covariance estimation 
library(xtable)             #for latex tables
library(boot)               #for bootstrapping
options(scipen = 999)       #forces R to use normal numbers instead of scientific notation


######################################################################
# Q1: Nonparametric bootstrap fail
######################################################################
set.seed(123)

N = 1000

# Simulate runif data
X = runif(N,0,1)

# Get max
x.max.obs = max(X)

# Write function for bootrap statistic
boot.stat = function(data, i){
  N*(x.max.obs-max(data[i]))
}

# Run bootsrap with 599 replications
boot.results = boot(data = X, R = 599, statistic = boot.stat)

# Make frequency plot
h         = hist(boot.results$t,plot=FALSE)
h$density = h$counts/sum(h$counts)
plot(h,freq=FALSE,main="Distribution of Bootstrap Statistic",xlab="Bootstrap statistic")


######################################################################
# Q2: Parametric bootstrap fail
######################################################################

# Generate parametric bootstrap samples
X.boot = replicate(599,runif(N,0,x.max.obs))

# Compute maximums for each replications
x.max.boot = sapply(1:599,function(i) max(X.boot[,i]))

# Compute bootstrap statistic
t.boot     = N*(x.max.obs-x.max.boot)

# Make frequency plot
h2         = hist(t.boot,plot=FALSE)
h2$density = h2$counts/sum(h2$counts)
plot(h2,freq=FALSE,main="Distribution of Parametric Bootstrap Statistic",xlab="Parametric bootstrap statistic",ylim=c(0,0.4),xlim=c(0,8))




