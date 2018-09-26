## ECON675: ASSIGNMENT 1
## Q3: ANALYSIS OF EXPERIMENTS
## Anirudh Yadav 
## 8/16/2018

######################################################################
# Load packages, clear workspace
######################################################################
rm(list = ls())             #clear workspace
library(dplyr)              #for data manipulation
library(ggplot2)            #for pretty plots
library(boot)               #for bootstrapping
options(scipen = 999)       #forces R to use normal numbers instead of scientific notation


######################################################################
# Input data
######################################################################

# Get LaLonde data
data <- read.csv('PhD_Coursework/ECON675/HW1/LaLonde_1986.csv')

# Convert to data frame
data <- as.data.frame(data)

######################################################################
# Q1 (a): difference in means estimator
######################################################################

# Rename variables
Y.obs       = data$earn78
treat       = data$treat

# Compute difference in means estimator
T.obs.dm    = mean(Y.obs[treat==1],na.rm=TRUE)-mean(Y.obs[treat==0],na.rm=TRUE)


######################################################################
# Q1 (b): conservative confidence intervals
######################################################################
N1 = sum(treat)
N0 = nrow(data)-N1

# Compute "conservative" standard error
s.1                <- sd(Y.obs[treat==1],na.rm=TRUE)^2
s.0                <- sd(Y.obs[treat==0],na.rm=TRUE)^2
se.conserv         <- sqrt(1/N1*s.1 + 1/N0*s.0)

# Compute lower and upper bounds of the interval and store in vector
CI.lower = T.obs.dm - qnorm(0.975)*se.conserv
CI.upper = T.obs.dm + qnorm(0.975)*se.conserv

# Store results
results            <- cbind(T.obs.dm,se.conserv,CI.lower,CI.upper)


######################################################################
# Q2 (a): Fisher Exact P-values
######################################################################

# The FEP function computes Fisher (approximate) p-values for 
# sharp null of no treatment effect, for the difference in means statistic
# and the K-S statistic.

# NOTES: 
# This function takes ~150 secs to run using the KS statistic with 250k draws!
# Is there a more efficient way to do this?


FEP <- function(K=249999,ks=FALSE){
  
  # Initialze vector of length K
  T.vec = vector(length=K)
  
  # Generate K random draws of the assignment vector
  T.MAT = replicate(K,sample(treat))
  
  if(!ks){
  
  # Compute observed difference in means
  T.obs    = mean(Y.obs[treat==1],na.rm=TRUE)-mean(Y.obs[treat==0],na.rm=TRUE)
  
      # Loop through random draws of the assignment vector, 
      # compute and store the test statistic
      for (i in 1:K) {
              T.dm    = mean(Y.obs[T.MAT[,i]==1],na.rm=TRUE)-mean(Y.obs[T.MAT[,i]==0],na.rm=TRUE)
              T.vec[i]            <- T.dm
            }
  
      }else{
    
    # USE K-S statistic
    options(warn=-1) #turn warnings off
    
    # Compute observed KS statistic
    T.obs               <- ks.test(Y.obs[treat==1],Y.obs[treat==0])$statistic
    
    # Loop through random draws of the assignment vector, 
    # compute and store the test statistic
    for (i in 1:K) {
      T.ks               <- ks.test(Y.obs[T.MAT[,i]==1],Y.obs[T.MAT[,i]==0])$statistic
      T.vec[i]           <- T.ks
    }
    
  }
  
  options(warn=0) #turn warnings back on!
  
  
  # Calculate p-value
  p = 1/K*sum(T.vec>=T.obs)
  
  return(p)
  
}


######################################################################
# Q2 (a): Fisher confidence intervals
######################################################################

## First I follow the approach in Imbens & Rubin, s5.7 ##

FisherInterval <- function(K=9999,C.vec=seq(5000,-1500,-250)){

  # Initialize vector of lenght C.vec
  P.vec = vector(length=length(C.vec))

  # Initialze vector of length K
  T.vec = vector(length=K)

  # Generate K random draws of the assignment vector
  T.MAT = replicate(K,sample(treat))
  
  for (j in 1:length(C.vec)){
    
    # Compute observed difference in means
    T.obs    = abs(mean(Y.obs[treat==1],na.rm=TRUE)-mean(Y.obs[treat==0],na.rm=TRUE)- C.vec[j])
    
    # Compute missing potential outcomes under the null
    Y.1 = ifelse(treat==1,Y.obs,Y.obs+C.vec[j])
    Y.0 = ifelse(treat==1,Y.obs-C.vec[j],Y.obs)

      for (i in 1:K) {
          T.dm    = abs(mean(Y.1[T.MAT[,i]==1],na.rm=TRUE)-mean(Y.0[T.MAT[,i]==0],na.rm=TRUE) - C.vec[j])
          T.vec[i]            <- T.dm
          }
  
    p = 1/K*sum(T.vec>=T.obs)
    P.vec[j] <- p
  }
  return(cbind(C.vec,P.vec))
}

# Run function with 10000 draws
# FisherInterval()


## Another way to compute the CI is using bootstrap ##

  # Compute missing potential outcomes under the null
  Y.1 = ifelse(treat==1,Y.obs,Y.obs+T.obs.dm)
  Y.0 = ifelse(treat==1,Y.obs-T.obs.dm,Y.obs)

  # Specify the statistic that we will compute for different permutations
  T.dm <- function(x, ind) {
      T.k <- mean(Y.1[data$treat[ind]==1]) - mean(Y.0[data$treat[ind]==0])
      return(T.k)
  }
  
  # Run bootstrap
  boot.result  <- boot(data = data, R = 9999, statistic = T.dm, sim = "permutation", stype = "i")
  boot.CI      <- quantile(boot.result$t, c(0.025, 0.975))
  
  # Empirical 95% CI for constant treatment effect = T.obs.dm
  print (boot.CI)

######################################################################
# Q3 (a): Plot power function
######################################################################

PowerFun <- function(x) {
            1 - pnorm(qnorm(0.975)-x/se.conserv) + pnorm(-qnorm(0.975)-x/se.conserv)
            }

# Plot usinging ggplot2
p1   <- ggplot(data.frame(x = c(-5000, 5000)), aes(x = x)) + stat_function(fun = PowerFun)

# Plot using curve
curve(1 - pnorm(qnorm(0.975)-x/se.conserv) + 
                pnorm(-qnorm(0.975)-x/se.conserv),-5000,5000,xlab="tau",ylab="Power")


######################################################################
# Q3 (b): Sample size calculation
######################################################################

# Parameterize the equation
p     = 2/3
tau   = 1000

# Write down the power function, which implicitly defines N
# [Note that I use the sample variances to proxy for the population variances]

Fun <- function(N){
        -0.8 + 1 - pnorm(qnorm(0.975)-tau/sqrt(1/N*s.1*(1/p)+1/N*s.0*(1/(1-p)))) +
          pnorm(-qnorm(0.975)-tau/sqrt(1/N*s.1*(1/p)+1/N*s.0*(1/(1-p)))) 
        }

# Solve for N
N.sol <- uniroot(Fun,c(0,100000000))$root

