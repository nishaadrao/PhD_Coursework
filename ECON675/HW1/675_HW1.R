## ECON675: ASSIGNMENT 1
## Q2: IMPLEMENTING LEAST SQUARES ESTIMATORS
## Anirudh Yadav 
## 8/16/2018

######################################################################
# Load packages, clear workspace
######################################################################
rm(list = ls())             #clear workspace
library('MASS')             #for ginv function
library('sandwich')         #for variance-covargiance estimation  
options(scipen = 999)       #forces R to use normal numbers instead of scientific notation


######################################################################
# Generate some random data
######################################################################
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


######################################################################
# Q4: Write function for computing OLS results + related stats
######################################################################

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
results2 <- linear_reg(X,Y,cholinv=TRUE)
diff     <- results1-results2

######################################################################
# Q5: Input data, create additional covariates
######################################################################

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

######################################################################
# Q5(a): run matrix implementation of OLS
######################################################################

# Run linear_reg function using different inverses
results1 <- linear_reg(X,Y)
results2 <- linear_reg(X,Y,cholinv=TRUE)


######################################################################
# Q5(b): compute OLS results using in-build lm() function
######################################################################

# Fit the linear regression
ols    <- lm(Y ~ X-1)      #nb. the minus 1 is becuase the lm package includes an intercept automatically
                           #    but the X matrix already includes an intercept

# Compute Eiker-White variance-covariance matrix (using sandwhich pkg)
# These match the standard errors from the manual computation above.
# Note that STATA's default when "r" is specified is equivalent to HC1!
V_check     <- vcovHC(ols, type = "HC0") 
se_check    <- as.matrix(sqrt(diag(V_check)))



