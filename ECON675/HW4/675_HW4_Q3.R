## ECON675: ASSIGNMENT 4
## Q3: POST-MODEL SELECTION INFERENCE
## Anirudh Yadav 
## 11/09/2018

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
options(scipen = 999)       #forces R to use normal numbers instead of scientific notation

######################################################################
# Generate random data
######################################################################

N     = 50
M     = 1000
SIGMA = matrix(c(1,0.85,0.85,1),2,2)

set.seed(1234)

# Generate covariates
W     = replicate(M,rmvnorm(N, mean = c(0,0), sigma = SIGMA, method="chol"))

# Generate errors
E     = replicate(M,rnorm(50))

# Generate outcomes
Y     = sapply(1:M,function(i) rep(1,N)+W[,,i]%*%c(0.5,1)+E[,i])

# Get beta.hats
beta.hats = sapply(1:M,function(i) lm(Y[,i]~W[,,i])$coefficients[2])

# Get t-stats for gamma.hats
t.stats   = sapply(1:M,function(i) summary(lm(Y[,i]~W[,,i]))[["coefficients"]][, "t value"][3])

# Get beta.tildes
beta.tildes = sapply(1:M,function(i) lm(Y[,i]~W[,1,i])$coefficients[2])

# Construct betas if the model selection is used
beta.sel    = ifelse(t.stats>=1.96,beta.hats,beta.tildes)
  
######################################################################
# [1] Summary Statistics for the different betas
######################################################################

# Summary statistics
beta.sum = rbind(summary(beta.hats),summary(beta.tildes),summary(beta.sel))

# Kernel density plots
plot(density(beta.hats,kernel="e",bw="ucv",na.rm=TRUE),main="Empirical Densities of the Different Estimators")
lines(density(beta.tildes,kernel="e",bw="ucv",na.rm=TRUE))
lines(density(beta.sel,kernel="e",bw="ucv",na.rm=TRUE))

# Make kernenl desity plot
plot.dat = data.frame(beta = c(beta.hats,beta.tildes,beta.sel),Estimator=rep(c("hat", "tilde","sel"), each = M))

densplot = ggplot(plot.dat,aes(x=beta,fill=Estimator))+ 
             geom_density(alpha=0.5, kernel="e",bw="ucv")+
             ggtitle("Kernel Density Plots for the Different Estimators")+
             xlab("Point Estimator")+
             ylab("Density")+
             theme(plot.title = element_text(hjust = 0.5))+
             scale_fill_discrete( 
                    name="Estimator",
                    breaks=c("hat", "tilde", "sel"),
                    labels=c("(i)", "(ii)", "(iii)"))+
             theme(legend.justification = c(0.05, 0.98), legend.position = c(0.05, 0.98))



