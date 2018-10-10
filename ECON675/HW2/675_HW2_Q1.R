## ECON675: ASSIGNMENT 2
## Q1: KERNEL DENSITY ESTIMATION
## Anirudh Yadav 
## 10/7/2018

######################################################################
# Load packages, clear workspace
######################################################################
rm(list = ls())             #clear workspace
library(dplyr)              #for data manipulation
library(data.table)         #for data manipulation
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

# Function for EP kernel
K.ep    <- function(x){
      y <- .75 * (1-x^2) * (abs(x) <= 1)
      return(y)
}

# Function to compute true density value
f.true  <- function(x){
     y<-0.5*dnorm(x,-1.5,sqrt(1.5))+0.5*dnorm(x,1,1)
     return(y)
}

# Create vector of bandwidths
h.list = h_aimse*seq(0.5,1.5,0.1)


# Generate vector of random draws from the given Gaussian DGP
N          <- 1000
M          <- 1000

components <- sample(1:2,prob=c(0.5,0.5),size=n,replace=TRUE)
mu.vec     <- c(-1.5,1)
sd.vec     <- sqrt(c(1.5,1))

randx     <- rnorm(n=N,mean=mu.vec[components],sd=sd.vec[components])


# Function for computing imses for a given bandwidth and random sample
imse         <- function(x.rand=randx, h=h_aimse){
  
  # First compute vector of density estimates at each x_i
  y   = sapply(x.rand,function(x) 1/(1000*h)*sum(K.ep((x.rand-x)/h)))
  
  # Convert y to data.table for easy manipulation
  y   = as.data.table(y)
  
  # Add true density values
  y[, y.true := f.true(x.rand)]
  
  # Compute squared errors
  y[, sq_er.li := (y - y.true)^2]
  
  # Compute imse.li
  imse.li <- y[, mean(sq_er.li)]
  
  # Compute leave-one-out fhats for each x_i
  lo   = sapply(1:N,function(i) 1/(1000*h)*sum(K.ep((as.matrix(x.rand)[-i,]-x.rand[i])/h)))
  
  # Add to data.table
  y[, fhat.lo := lo]
  
  # Compute squared errors
  y[, sq_er.lo := (fhat.lo - y.true)^2]
  
  # Compute imse
  imse.li <- y[, mean(sq_er.li)]
  
  # Compute imse.lo
  imse.lo <- y[, mean(sq_er.lo)]
  
  output <- c(imse.li,imse.lo)
  
  return(output)
}  

# Compute mse for each h in h.list
imse.vec      <- sapply(h.list,imse)

# THIS IS CLOSE TO BEING VERY EFFICIENT!!!
# JUST NEED TO FIGURE OUT HOW TO LOOP OVER THE h's and the X's together!!
X.mat      <- replicate(M,rnorm(n=N,mean=mu.vec[components],sd=sd.vec[components]))
system.time(yo<-sapply(1:M,function(i) imse(X.mat[,i])))


  