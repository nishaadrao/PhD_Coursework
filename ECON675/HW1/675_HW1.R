## ECON675: ASSIGNMENT 1
## Q2: IMPLEMENTING LEAST SQUARES ESTIMATORS
## Anirudh Yadav 
## 8/16/2018

#-------Load packages-----------------------------------------
#--------------------------------------------------------------
rm(list = ls())             #clear workspace
library('MASS')             #for ginv function
library('sandwich')         #for variance-covargiance estimation  
options(scipen = 999)       #forces R to use normal numbers instead of scientific notation

#-------Input data & data processing --------------------------
#--------------------------------------------------------------

# Get LaLonde data
data <- read.csv('PhD_Coursework/ECON675/HW1/LaLonde_1986.csv')

# Convert to data frame
df <- as.data.frame(data)

# Add a column of ones for intercept in regression
df$ones <- 1

# Create squared education variable
df$educsq <- (df$educ)^2

# Create black-earn74 interaction 
df$black_earn74 <- (df$black)*(df$earn74)


#-------OLS point estimator------------------------------------
#--------------------------------------------------------------

# Create X matrix
X <- cbind(df$ones,df$treat,df$black,df$age,df$educ,df$educsq,df$earn74,df$black_earn74,df$u74,df$u75)

# Create Y vector
Y <- cbind(df$earn78)

# Compute OLS point estimator
beta <- solve(crossprod(X))%*%crossprod(X,Y)
rownames(beta) <- c("const","treat","black","age","educ","educsq","earn74","black_earn74","u74","u75")


#-------Variance-covariance estimator--------------------------
#--------------------------------------------------------------

# Compute residuals vector
e     <- Y - X%*%beta

# Construct diagonal matrix of squared residuals 
D     <- diag(as.numeric(e^2))

# Compute Eiker-White variance-covariance matrix (finite sample version)
# NOTE: this corresponds to eqn 4.37 in Hansen (p 99), w/o df adjustment
n <- nrow(df)
k <- ncol(X)

V     <- solve(crossprod(X))%*%t(X)%*%D%*%X%*%solve(crossprod(X))

# Compute standard errors
se    <- as.matrix(sqrt(diag(V)))

#------- t-stats and p-values----------------------------------
#--------------------------------------------------------------

# Compute t-stats for testing H0:beta_k=0 
t=matrix()
for (i in 1:10) {
  t[i] = beta[i]/se[i]
  }
t <- as.matrix(t)

# Compute p-values
p=matrix()
for (i in 1:10) {
  p[i] = 2*pt(abs(t[i]),df=n-k,lower.tail=FALSE)
  }
p <- as.matrix(p)


#------- 95 pct confidence intervals---------------------------
#--------------------------------------------------------------

#Compute lower bounds of the CI
CIlower=matrix()
for (i in 1:10) {
  CIlower[i] = beta[i] - qnorm(0.975)*se[i]
  }

#Compute upper bounds of the CI
CIupper=matrix()
for (i in 1:10) {
  CIupper[i] = beta[i] + qnorm(0.975)*se[i]
  }


#------- Combine results into a dataframe----------------------
#--------------------------------------------------------------
results             <- cbind(beta,se,t,p,as.matrix(CIlower),as.matrix(CIupper))
colnames(results)   <- c("beta","se","t","p","CIlower","CIupper")

#-------OLS with lm function-----------------------------------
#--------------------------------------------------------------
# Below I use R's lm() function to run the regression above (Q5(b))

# Fit the linear regression
ols    <- lm(Y ~ X-1)      #nb. the minus 1 is becuase the lm package includes an intercept automatically
                           #    but the X matrix already includes an intercept

# Compute Eiker-White variance-covariance matrix (using sandwhich pkg)
V_check     <- vcovHC(ols, type = "HC") 

# Compute SEs
se_check    <- as.matrix(sqrt(diag(V_check)))


