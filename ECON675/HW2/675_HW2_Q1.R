## ECON675: ASSIGNMENT 2
## Q1: KERNEL DENSITY ESTIMATION
## Anirudh Yadav 
## 10/7/2018

######################################################################
# Load packages, clear workspace
######################################################################
rm(list = ls())             #clear workspace
library(dplyr)              #for data manipulation
library(ggplot2)            #for pretty plots
library(boot)               #for bootstrapping
options(scipen = 999)       #forces R to use normal numbers instead of scientific notation



######################################################################
# Q3 (a): compute theoretically optimal bandwidth
######################################################################
# NB. This code only makes sense with the associated tex file...

# Write function to compute second derivative of normal density
d2norm  <- function(x, mu=0, v=1) {
  dnorm(x,mu,sqrt(v))*(((x-mu)/v)^2-1/v)
}

# Second derivative, squared of given Gaussian mixture
myf     <- function(x){
  (0.5*d2norm(x,-1.5,1.5)+0.5*d2norm(x,1,1))^2
}

# Compute required integral
k1     <-integrate(myf, -Inf, Inf)$val

# Compute optimal bandwidth
n      <- 1000
k2     <- 1/5
k3     <- 3/5
P      <- 2

h_aimse <- ((1/(2*P*n))*(factorial(P)/k2)^2*(k3/k1))^(1/(1+2*P))


######################################################################
# Q3 (b): monte carlo
######################################################################
set.seed(0)

# Write function for EP kernel
K <- function(x){
  y <- .75 * (1-x^2) * (abs(x) <= 1)
}

# Generate equally spaced points at which the density is to be estimated
x.grid     <- seq.int(from=-7, to=7, length.out = 1000)

# Generate big matrix of random draws from the given Gaussian DGP
N          <- 1000
M          <- 1000

components <- sample(1:2,prob=c(0.5,0.5),size=n,replace=TRUE)
mu.vec     <- c(-1.5,1)
sd.vec     <- sqrt(c(1.5,1))

x.rand     <- rnorm(n=N,mean=mu.vec[components],sd=sd.vec[components])

# Each column of X.mat contains 1000 draws from the DGP
X.mat      <- replicate(M,rnorm(n=N,mean=mu.vec[components],sd=sd.vec[components]))

# Create vector of bandwidths
h.list = h_aimse*seq(0.5,1.5,0.1)

# Compute density estimates for a given bandwidth
fhat       <- function(h=h_aimse){
  sapply(x.grid,function(x) 1/(1000*h)*sum(K((x.rand-x)/h)))
}  

fhatmat       <- function(h=h_aimse){
  sapply(x.grid,function(x) 1/(1000*h)*colSums(K((X.mat-x)/h)))
} 

# Compute density estimates for each h in h.list
y.mat      <- sapply(h.list,fhat)

#y.mat.big<- sapply(h.list,fhatmat)
  