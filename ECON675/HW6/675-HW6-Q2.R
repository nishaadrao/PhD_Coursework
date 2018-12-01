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
library(rdrobust)           #for RD plots and other stuff
library(rddensity)          #for RD density continuity tests
options(scipen = 999)       #forces R to use normal numbers instead of scientific notation

######################################################################
# Input data
######################################################################
setwd("/Users/Anirudh/Desktop/GitHub/PhD_Coursework/ECON675/HW6")
data <- as.data.table(read.csv('HeadStart.csv'))

######################################################################
# [2.1] RD Plots of pre-intervention mortality
######################################################################

# Evenly-spaced bins, IMSE optimal
rdplot(data[,mort_related_pre],data[,povrate60],p=1,binselect = "es",x.label="povrate60",y.label="mort_related_pre",title="")
dev.copy(pdf,'q2-1-es.pdf')
dev.off()

# Evenly-spaced bins, mimicking variance 
rdplot(data[,mort_related_pre],data[,povrate60],p=1,binselect = "esmv",x.label="povrate60",y.label="mort_related_pre",title="")
dev.copy(pdf,'q2-1-esmv.pdf')
dev.off()

# Quantile-spaced bins, IMSE optimal
rdplot(data[,mort_related_pre],data[,povrate60],p=1,binselect = "qs",x.label="povrate60",y.label="mort_related_pre",title="")
dev.copy(pdf,'q2-1-qs.pdf')
dev.off()

# Quantile-spaced bins, mimicking variance
rdplot(data[,mort_related_pre],data[,povrate60],p=1,binselect = "qsmv",x.label="povrate60",y.label="mort_related_pre",title="")
dev.copy(pdf,'q2-1-qsmv.pdf')
dev.off()

######################################################################
# [2.1] Formal falsification tests
######################################################################

## Exact binomial tests for different windows around the cutoff

# Vector of windows
h.vec = seq(0.3,1.3,0.2)

# Get running variable
x     = data[,povrate60]

# Number of observations just above and below the cutoff
N.l   = sapply(1:length(h.vec),function(i) sum(x >= -h.vec[i] & x <=0))
N.u   = sapply(1:length(h.vec),function(i) sum(x >= 0 & x <= h.vec[i]))

# Total number of observations in the window
N.t   = N.l + N.u

# Conduct exact binomial tests (p=0.5), where success is treatment and store p-vals
binom.pvals = sapply(1:length(h.vec),function(i) binom.test(N.u[i],N.t[i])$p.value)

# Put results together for latex
binom.results = cbind(h.vec,N.l,N.u,binom.pvals)
xtable(binom.results,digits = c(0,1,0,0,3))

## Continuity in density tests
