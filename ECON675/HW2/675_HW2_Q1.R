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
# NB. This code only makes sense with the associated tex file

# Write function to compute second derivative of normal density
d2norm  <- function(x, mu=0, v=1) {
  dnorm(x)*(((x-mu)/v)^2-v)
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

h      <- ((1/2*P*n)*(factorial(P)/k2)^2*(k3/k1))^(1/(1+2*P))






