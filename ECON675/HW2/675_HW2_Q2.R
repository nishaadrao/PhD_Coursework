## ECON675: ASSIGNMENT 2
## Q2: LINEAR SMOOTHERS, CROSS-VALIDATION AND SERIES
## Anirudh Yadav 
## 10/10/2018

######################################################################
# Load packages, clear workspace
######################################################################
rm(list = ls())             #clear workspace
library(foreach)            #for looping
library(dplyr)              #for data manipulation
library(data.table)         #for data manipulation
library(ggplot2)            #for pretty plots
library(boot)               #for bootstrapping
library(Matrix)             #fast matrix calcs
options(scipen = 999)       #forces R to use normal numbers instead of scientific notation


######################################################################
# Q5 (a): generate random data
######################################################################
N   = 1000
M   = 1000

# Generate x's
x.mat   = replicate(M,runif(N,-1,1))

# Generate chi-squared r.v
chi = replicate(M,rchisq(N,5))

# Generate errors
u.mat   = (chi-5)*x.mat

# Generate y's
y.mat   = exp(-0.1*(4*x.mat-1)^2)*sin(5*x.mat) + u.mat

######################################################################
# Q5 (b): cross-validation series estimator
######################################################################


cross.val <- function(i){

data = data.table(y=y.mat[,i],x=x.mat[,i],const=1)
temp = rep(NaN, 20)

for (k in 1:20) {
  
  data[, temp := x^k]
  setnames(data, "temp", paste0("x_", k))
  
  X <- as.matrix(data[, c("const",grep("x_", colnames(data), value = TRUE)), with = FALSE])
  Y <- as.matrix(data[,y])
  
  X.Q   <- qr.Q(qr(X))
  XX <- X.Q %*% t(X.Q)
  Y.hat <- XX %*% Y
  W <- diag(XX)
  
  temp[k] <- mean(((Y-Y.hat) / (1-W))^2)
  
}
return(temp)
}

# RUN SIMULATION -- RUNTIME 5 MINS
results <- sapply(1:M,function(i) cross.val(i))

results.avg=rowMeans(results)

K.hat=which.min(results.avg)


# #Plot CV
# g <- as.data.frame(cbind(1:20,results.avg))
# colnames(g) <- c("K", "CV")
# 
# 
# ggplot(g,aes(x=K, y=CV)) +
#   geom_line(linetype = "dashed")+
#   geom_point()+
#   labs(title="Simulated Cross-Validation Errors for M=1000 Simulations") +theme(plot.title = element_text(hjust = 0.5))

######################################################################
# Q5 (c): diagnostics
######################################################################

# Generate grid of x-valus
x.vec = seq(-1,1,0.1)

# Write function to compute optimal beta's (i.e. K=7) and standard errors
cv.beta <- function(i){
  
  data = data.table(y=y.mat[,i],x=x.mat[,i],const=1)
  
  for (k in 1:7) {
    
    data[, temp := x^k]
    setnames(data, "temp", paste0("x_", k))
  }  
    
    X <- as.matrix(data[, c("const",grep("x_", colnames(data), value = TRUE)), with = FALSE])
    Y <- as.matrix(data[,y])
    
    beta <- solve(crossprod(X))%*%crossprod(X,Y) 
    
    X.Q   <- qr.Q(qr(X))
    XX <- X.Q %*% t(X.Q)
    Y.hat <- XX %*% Y
    W <- diag(XX)
    
    return(t(beta)) 

}

# Write function to compute optimal standard errors
cv.se <- function(i){
  
  data = data.table(y=y.mat[,i],x=x.mat[,i],const=1)
  
  for (k in 1:7) {
    
    data[, temp := x^k]
    setnames(data, "temp", paste0("x_", k))
  }  
  
  X <- as.matrix(data[, c("const",grep("x_", colnames(data), value = TRUE)), with = FALSE])
  Y <- as.matrix(data[,y])
  
  X.Q   <- qr.Q(qr(X))
  XX <- X.Q %*% t(X.Q)
  Y.hat <- XX %*% Y
  W <- diag(XX)
  
  s <- (1/(N-1))*sum((Y-Y.hat)^2)
  
  V <- s*(t(W)%*%W)
  
  se <- sqrt(V)
  
  return(se) 
  
}

# Get optimal betas for each m
results.beta <- sapply(1:M,function(i) cv.beta(i))

# Get associated standard errors
results.se <- sapply(1:M,function(i) cv.se(i))

# Compute averages over M
opt.beta  <- rowMeans(results.beta)
opt.se    <- rowMeans(results.se)

# Compute regressors for each number in the x.grid
X.new     <- t(sapply(x.vec, function(x) return(cbind(1,x,x^2,x^3,x^4,x^5,x^6,x^7))))

# Compute y.hats
y.hats    <- X.new%*%as.vector(opt.beta)

# Write the true regression function
f.true <- function(x) exp(-0.1*(4*x-1)^2)*sin(5*x)

# MAKE PLOT

plot.data = as.data.frame(cbind(x.vec,y.hats))

ggplot(plot.data,aes(x=x.vec,y=V2))+
  geom_line(linetype = "dashed")+geom_point()+
  stat_function(fun =f.true)+
  labs(title="True and Series Estimate of the Regression function",y="y",x="x") +theme(plot.title = element_text(hjust = 0.5))






# 