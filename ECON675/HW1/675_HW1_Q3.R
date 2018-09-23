## ECON675: ASSIGNMENT 1
## Q3: ANALYSIS OF EXPERIMENTS
## Anirudh Yadav 
## 8/16/2018

#------- Load packages ----------------------------------------#
#______________________________________________________________#
rm(list = ls())             #clear workspace
library('dplyr')            #for data manipulation
options(scipen = 999)       #forces R to use normal numbers instead of scientific notation


#------- Input data & processing ------------------------------#
#______________________________________________________________#

# Get LaLonde data
data <- read.csv('PhD_Coursework/ECON675/HW1/LaLonde_1986.csv')

# Convert to data frame
data <- as.data.frame(data)

# Attach data frame so we can refer to variables just using their names
attach(data)

#==============================================================#
#================= QUESTION 1: NEYMAN's APPROACH ==============#
#==============================================================#
Y.obs       = earn78

#------- Compute difference in means estimator ----------------#
#______________________________________________________________#
T.obs.dm    = mean(Y.obs[treat==1],na.rm=TRUE)-mean(Y.obs[treat==0],na.rm=TRUE)


#------- Construct conservative CIs  --------------------------#
#______________________________________________________________#
N1 = sum(treat)
N0 = nrow(data)-N1

# Compute "conservative" standard error
se.conserv         <- sqrt(1/N1*sd(Y.obs[treat==1],na.rm=TRUE)^2 + 1/N0*sd(Y.obs[treat==0],na.rm=TRUE)^2)

# Compute lower and upper bounds of the interval and store in vector
CI.lower = T.obs.dm - qnorm(0.975)*se.conserv
CI.upper = T.obs.dm + qnorm(0.975)*se.conserv

results            <- cbind(T.obs.dm,se.conserv,CI.lower,CI.upper)


#==============================================================#
#================= QUESTION 2: FISHER'S APPROACH ==============#
#==============================================================#

#------- P-values ---------------------------------------------#
#______________________________________________________________#


# The FEP function computes Fisher (approximate) p-values for 
# sharp null of no treatment effect, for the difference in means statistic
# and the K-S statistic.

# NOTES: 
# [1] This function takes ~150 secs to run using the KS statistic with 250k draws!
# Is there a more efficient way to do this?
# [2] Below, I've used the in-built ks.test() function to get the K-S statistic; 
# I tried to look at the source code to see how to compute the statistic manually, but it was
# too hard!

FEP <- function(K=250000,ks=FALSE){
  
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

#------- 95% Confidence Interval-------------------------------#
#______________________________________________________________#


FisherInterval <- function(K=10000,C.vec=seq(10000,-1500,-250),ks=FALSE){

  # Initialize vector of lenght C.vec
  P.vec = vector(length=length(C.vec))

  # Initialze vector of length K
  T.vec = vector(length=K)

  # Generate K random draws of the assignment vector
  T.MAT = replicate(K,sample(treat))

  # Compute observed difference in means
  T.obs    = mean(Y.obs[treat==1],na.rm=TRUE)-mean(Y.obs[treat==0],na.rm=TRUE)
  
  for (j in 1:length(C.vec)){

    # Compute missing potential outcomes under the null
    Y.1 = ifelse(treat==1,earn78,earn78+C.vec[j])
    Y.0 = ifelse(treat==1,earn78-C.vec[j],earn78)

      for (i in 1:K) {
          T.dm    = mean(Y.1[T.MAT[,i]==1],na.rm=TRUE)-mean(Y.0[T.MAT[,i]==0],na.rm=TRUE) - C.vec[j]
          T.vec[i]            <- T.dm
          }
  
    p = 1/K*sum(T.vec>=T.obs)
    P.vec[j] <- p
  }
  return(cbind(C.vec,P.vec))
}


