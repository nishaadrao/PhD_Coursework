## ECON675: ASSIGNMENT 1
## Q2: IMPLEMENTING LEAST SQUARES ESTIMATORS
## Q3: ANALYSIS OF EXPERIMENTS
## Anirudh Yadav 
## 8/16/2018

#-------Load packages-----------------------------------------
#-------------------------------------------------------------
rm(list = ls())             #clear workspace
library('MASS')             #for ginv function
library('sandwich')         #for variance-covargiance estimation  
options(scipen = 999)       #forces R to use normal numbers instead of scientific notation


#==============================================================#
#================= QUESTION 2: OLS ============================#
#==============================================================#


##################### PART 4 ###################################


#------- Generate some random data ----------------------------#
#______________________________________________________________#

# Draw some random numbers from the std norm distribution
# and use these as our data for testing the OLS function below.
# Specifically, we'll use these data to test whether the two
# different types of inverses give the same results.

# Generate covariates
X        <- replicate(3, rnorm(100))

# Generate iid mean-zero errors
U        <- rnorm(100)

# Impose a `true' beta
betatrue <- as.matrix(c(1,2,3))

# Compute the resulting Y
Y        <- X%*%betatrue + U


#------- Write function for OLS results -----------------------#
#______________________________________________________________#

linear_reg <- function(X,Y, cholinv=FALSE, alpha=0.05){

# Compute crossproduct matrix 
  
  if(!cholinv){
      # using symmetric inverse
      M = solve(crossprod(X))
      
  }else{
      #using Cholesky inverse
      m  = chol(crossprod(X))
      M = chol2inv(m)
  }
    
      # Compute OLS estimator
      beta = M%*%crossprod(X,Y)  

      # Compute residuals vector
      u     <- Y - X%*%beta

      # Construct diagonal matrix of squared residuals 
      D     <- diag(as.numeric(u^2))

      # Compute variance-covarinance matrix
      n <- nrow(X)
      k <- ncol(X)
      V     <- n*M%*%t(X)%*%D%*%X%*%M

      # Compute standard errors
      se    <- as.matrix(sqrt(diag(V)/n))

      # Compute t-stats for testing H0:beta_k=0 
      t=matrix()
      for (i in 1:k) {
              t[i] = beta[i]/se[i]
      }
      t <- as.matrix(t)

      # Compute p-values
      p=matrix()
      for (i in 1:k) {
              p[i] = 2*pt(abs(t[i]),df=n-k,lower.tail=FALSE)
      }
      p <- as.matrix(p)

      #Compute lower bounds of the CI
      CIlower=matrix()
      for (i in 1:k) {
               CIlower[i] = beta[i] - qnorm(1-alpha/2)*se[i]
      }

      #Compute upper bounds of the CI
      CIupper=matrix()
      for (i in 1:k) {
              CIupper[i] = beta[i] + qnorm(1-alpha/2)*se[i]
      }
  
  
    #Combine results into a dataframe
    results             <- cbind(beta,se,t,p,as.matrix(CIlower),as.matrix(CIupper))
    colnames(results)   <- c("beta","se","t","p","CIlower","CIupper")
  
  return(results)
}

# Test whether choice of inverse affects results
results1 <- linear_reg(X,Y)
results2 <- linear_reg(X,Y)

diff <- results1-results2

##################### PART 5 ###################################


#------- Input data & processing ------------------------------#
#______________________________________________________________#

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

# Create X matrix
X <- cbind(df$ones,df$treat,df$black,df$age,df$educ,df$educsq,df$earn74,df$black_earn74,df$u74,df$u75)

# Create Y vector
Y <- cbind(df$earn78)


#-------  Run linear_reg function -----------------------------#
#______________________________________________________________#

# Run linear_reg function using different inverses
results1 <- linear_reg(X,Y)
results2 <- linear_reg(X,Y,cholinv=TRUE)



#-------  Compute OLS results with lm() -----------------------#
#______________________________________________________________#

# Below I use R's lm() function to run the regression above (Q5(b))

# Fit the linear regression
ols    <- lm(Y ~ X-1)      #nb. the minus 1 is becuase the lm package includes an intercept automatically
                           #    but the X matrix already includes an intercept

# Compute Eiker-White variance-covariance matrix (using sandwhich pkg)
# These match the standard errors from the manual computation above.
# Note that STATA's default when "r" is specified is equivalent to HC1!
V_check     <- vcovHC(ols, type = "HC0") 
se_check    <- as.matrix(sqrt(diag(V_check)))



#==============================================================#
#================= QUESTION 3: EXPERIMENTS ====================#
#==============================================================#


################### PART 1: NEYMAN'S APPROACH ##################

# (a): Compute difference in means estimator
N1 = sum(df$treat)
N0 = nrow(df)-N1

meanearnings       <- as.vector(tapply(df$earn78,df$treat,mean))
Tdm                <- meanearnings[2]-meanearnings[1]

# (b): Construct conservative 95% CI
N1 = sum(df$treat)
N0 = nrow(df)-N1
  
  # Get sample standard deviations of earn78 by treatment group
  stddevs          <- aggregate(df$earn78,list(treat=df$treat),sd)
  # write some manual code to check that this is using the correct sd formula
  samplevars       <- as.matrix(stddevs^2)
  
  seconserv         <- sqrt(1/N1*samplevars[2,2] + 1/N0*samplevars[1,2])
  
  # Compute lower and upper bounds of the interval and store in vector
  CIlower = Tdm - qnorm(0.975)*seconserv
  CIupper = Tdm + qnorm(0.975)*seconserv
  
results3            <- cbind(Tdm,seconserv,CIlower,CIupper)
  