###################################################################
# ECON 675, Assignment 2
# Fall 2018
###################################################################

###################################################################
# Q1. 3. (a)
###################################################################

# parameters
mu1     <- -1.5
mu2     <- 1
sigma1  <- sqrt(1.5)
sigma2  <- sqrt(1)
prob    <- 0.5
n       <- 1000
M       <- 100

# The Hermite polynomial (for computing the derivatives)
hermite <- function(x, p){
  if (p==0)  out = 1
  if (p==1)  out = x
  if (p==2)  out = x^2 - 1
  if (p==3)  out = x^3 - 3*x
  if (p==4)  out = x^4 - 6*x^2 + 3
  if (p==5)  out = x^5 - 10*x^3 + 15*x
  if (p==6)  out = x^6 - 15*x^4 + 45*x^2 - 15
  if (p==7)  out = x^7 - 21*x^5 + 105*x^3 - 105*x
  if (p==8)  out = x^8 - 28*x^6 + 210*x^4 - 420*x^2 + 105
  if (p==9)  out = x^9 - 36*x^7 + 378*x^5 - 1260*x^3 + 945*x
  if (p==10) out = x^10 - 45*x^8 + 630*x^6 - 3150*x^4 + 4725*x^2 - 945
  return(out)
}

# The kernel function
K_0 <- function(x) {
  return( 0.75 * (1-x^2) * (abs(x)<=1) )
}

# derivative of kernel
K_1 <- function(x) {
  return( 0.75 * (-2 * x) * (abs(x)<=1) )
}

# The normal density and its derivatives
f_norm <- function(x, mu, sigma, s) {
  y <- (x-mu) / sigma
  temp <- dnorm(y) * hermite(y, s) / sigma^(s+1)
  return(temp)
}

# The DGP density and its derivatives
f_DGP <- function(x, s) {
  temp <- prob * f_norm(x, mu1, sigma1, s) +
    (1 - prob) * f_norm(x, mu2, sigma2, s)
  return(temp)
}

# The `mu' functional
mu_f <- function(f, l, low=-1, up=1) {
  integrand <- function(x) { x^l * f(x) }
  return( integrate(integrand, lower=low, upper=up)$value )
}

# the `theta' functional
theta_f <- function(f, low=-1, up=1) {
  integrand <- function(x) { f(x)^2 }
  return( integrate(integrand, lower=low, upper=up)$value )
}

# The asymptotic IMSE (AIMSE) 
AIMSE <- function(h, s) {
  if (s == 0) {
    K <- function(x) { return(K_0(x)) }
    f <- function(x) { return(f_DGP(x, s=2)) }
  } else {
    K <- function(x) { return(K_1(x)) }
    f <- function(x) { return(f_DGP(x, s=3)) }
  }
  variance <- theta_f(K) / (n*h^(2*s+1))
  bias <- theta_f(f, low=-10, up=10) * mu_f(K_0, 2)^2 * h^(4) / 4 
  return( variance + bias )
}

# The AIMSE-optimal bandwidth 
h_AIMSE <- function(s) {
  if (s == 0) {
    K <- function(x) { return(K_0(x)) }
    f <- function(x) { return(f_DGP(x, s=2)) }
  } else {
    K <- function(x) { return(K_1(x)) }
    f <- function(x) { return(f_DGP(x, s=3)) }
  }
  num <- theta_f(K) * (2*s + 1)
  denom <- theta_f(f, low=-10, up=10) * mu_f(K_0, 2)^2 * n
  return( (num/denom)^(1/(2*s+5)))
}

# Two ways to get the AIMSE-optimal bandwidth
h_AIMSE(0)
optimize(AIMSE, c(0.01, 2), s=0)$minimum
h_AIMSE(1)
optimize(AIMSE, c(0.01, 2), s=1)$minimum

###################################################################
# Q1. 3. (b)
###################################################################
# DGP
data_gen <- function() {
  temp <- rep(NaN, n)
  N1 <- rbinom(1, size=n, prob=prob)
  temp[1:N1] <- rnorm(N1, mean=mu1, sd=sigma1)
  temp[(N1+1):n] <- rnorm(n-N1, mean=mu2, sd=sigma2)
  return(temp)
}

h_AIMSE_opt <- h_AIMSE(0)
h <- h_AIMSE_opt * seq(from=0.5, to=1.5, by=0.1)
nh <- length(h)

# The estimator
f_hat <- function(x, X, h) {
  temp <- K_0((X-x) / h) / h
  return( mean(temp) )
}

IMSE_LI <- matrix(NaN, nrow=M, ncol=nh)
IMSE_LO <- IMSE_LI

set.seed(123)
ptm <- proc.time()
for (m in 1:M) {
  X <- data_gen()
  for (j in 1:nh) {
    temp_LI <- rep(NaN, n)
    temp_LO <- temp_LI
    for (i in 1:n) {
      temp_LI[i] <- (f_hat(X[i], X, h[j]) - f_DGP(X[i], 0))^2
      temp_LO[i] <- (f_hat(X[i], X[-i], h[j]) - f_DGP(X[i], 0))^2
    }
    IMSE_LI[m, j] <- mean(temp_LI)
    IMSE_LO[m, j] <- mean(temp_LO)
  }
}
proc.time() - ptm


plot(h, colMeans(IMSE_LI), type="l", xlab="bandwidth", ylab="IMSE", main="")
lines(h, colMeans(IMSE_LO), type="l", lty=2)
lines(c(h_AIMSE(0), h_AIMSE(0)), c(-1, 1), type="l", col="red", lwd=3)
lines(c(h_AIMSE(0), h[which.min(colMeans(IMSE_LI))]), c(-1, 1), type="l", lty=1)
lines(c(h_AIMSE(0), h[which.min(colMeans(IMSE_LO))]), c(-1, 1), type="l", lty=2)

h[which.min(colMeans(IMSE_LI))]
h[which.min(colMeans(IMSE_LO))]

###################################################################
# Q1. 3. (d)
###################################################################

h_hat <- function(X) {
  mu_hat <- mean(X)
  sigma_hat <- sd(X)
  f <- function(x) { return(f_norm(x, mu_hat, sigma_hat, 2)) }
  num <- theta_f(K_0) 
  denom <- theta_f(f, low=mu_hat-6*sigma_hat, up=mu_hat+6*sigma_hat) * mu_f(K_0, 2)^2 * n
  return( (num/denom)^(1/5) )
}

h_temp <- rep(NaN, M)

set.seed(123)
ptm <- proc.time()
for (m in 1:M) {
  X <- data_gen()
  h_temp[m] <- h_hat(X)
}
proc.time() - ptm

hist(h_temp, breaks=40, freq=FALSE, xlab="bandwidth", main="")
lines(c(mean(h_temp), mean(h_temp)), c(-1, 50), col="red", lwd=3)

summary(h_temp)

###################################################################
# Q3. 4. (b)
###################################################################

# generate the polynomial basis
gen.P = function(Z,K) {
  if (K==0)   out = NULL;
  if (K==1)   out = poly(Z,degree=1,raw=TRUE);
  if (K==2)  {out = poly(Z,degree=1,raw=TRUE); for (j in 1:ncol(Z)) out = cbind(out,Z[,j]^2);}
  if (K==2.5) out = poly(Z,degree=2,raw=TRUE);
  if (K==3)  {out = poly(Z,degree=2,raw=TRUE); for (j in 1:ncol(Z)) out = cbind(out,Z[,j]^3);}
  if (K==3.5) out = poly(Z,degree=3,raw=TRUE);
  if (K==4)  {out = poly(Z,degree=3,raw=TRUE); for (j in 1:ncol(Z)) out = cbind(out,Z[,j]^4);}
  if (K==4.5) out = poly(Z,degree=4,raw=TRUE);
  if (K==5)  {out = poly(Z,degree=4,raw=TRUE); for (j in 1:ncol(Z)) out = cbind(out,Z[,j]^5);}
  if (K==5.5) out = poly(Z,degree=5,raw=TRUE);
  if (K>=6)  {out = poly(Z,degree=5,raw=TRUE); for (k in 6:K) for (j in 1:ncol(Z)) out = cbind(out,Z[,j]^k);}
  ## RETURN POLYNOMIAL BASIS
  return(out)
}

# DGP
data_gen <- function(n) {
  X <- matrix((runif(n*5)-0.5)*2, ncol=5)
  V <- matrix(rnorm(n), ncol=1)
  G <- matrix(exp(diag(X %*% t(X))), ncol=1)
  E <- matrix(0.3637899 * (1 + diag(X %*% t(X))) * V, ncol=1)
  U <- matrix(rnorm(n), ncol=1)
  tt <- matrix(sqrt(diag(X %*% t(X))) + U >= 0, ncol=1) * 1
  Y <- matrix(tt + G + E, ncol=1)
  return(list(Y=Y, X=X, tt=tt))
} 

n   <- 500
K   <- c(1, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6, 7, 8, 9, 10)
K.r <- c(6, 11, 21, 26, 56, 61, 126, 131, 252, 257, 262, 267, 272, 277)
nK  <- length(K)
M   <- 1000
theta.hat <- matrix(NaN, ncol=nK, nrow=M)
se.hat    <- theta.hat

set.seed(123)
ptm <- proc.time()
for (m in 1:M) {
  data <- data_gen(n)
  X <- data$X; Y <- data$Y; tt <- data$tt
  for (k in 1:nK) {
    X.pol <- cbind(1, gen.P(X, K[k]))
    X.Q   <- qr.Q(qr(X.pol))
    MP     <- diag(rep(1,n)) - X.Q %*% t(X.Q)
    Y.M <- MP %*% Y
    tt.M <- MP %*% tt
    theta.hat[m, k] <- (t(tt.M) %*% Y.M) / (t(tt.M) %*% tt.M)
    Sigma <- diag((as.numeric((Y.M - tt.M*theta.hat[m, k])))^2)
    se.hat[m, k] <- sqrt(t(tt.M) %*% Sigma %*% tt.M) / (t(tt.M) %*% tt.M)
  }
}
proc.time() - ptm

table <- matrix(NaN, ncol=6, nrow=nK)
for (k in 1:nK) {
  table[k, 1] <- K.r[k]
  table[k, 2] <- mean(theta.hat[, k]) - 1                           # bias
  table[k, 3] <- sd(theta.hat[, k])                                 # standard deviation
  table[k, 4] <- table[k, 2]^2 + table[k, 3]^2                      # mse
  table[k, 5] <- mean(se.hat[, k])                                  # mean standard error
  table[k, 6] <- mean((theta.hat[, k] - 1.96 * se.hat[, k] > 1) | 
                        (theta.hat[, k] + 1.96 * se.hat[, k] < 1))  # rejection rate
}
write.table(round(table,3), "partial_linear.txt", sep=" & ", eol="\\\\ \n", col.names = FALSE, row.names = FALSE)

###################################################################
# Q3. 4. (c)
###################################################################

n   <- 500
K   <- c(1, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6, 7, 8, 9, 10)
K.r <- c(6, 11, 21, 26, 56, 61, 126, 131, 252, 257, 262, 267, 272, 277)
nK  <- length(K)
M   <- 1000

# generate the polynomial basis
gen.P = function(Z,K) {
  if (K==0)   out = NULL;
  if (K==1)   out = poly(Z,degree=1,raw=TRUE);
  if (K==2)  {out = poly(Z,degree=1,raw=TRUE); for (j in 1:ncol(Z)) out = cbind(out,Z[,j]^2);}
  if (K==2.5) out = poly(Z,degree=2,raw=TRUE);
  if (K==3)  {out = poly(Z,degree=2,raw=TRUE); for (j in 1:ncol(Z)) out = cbind(out,Z[,j]^3);}
  if (K==3.5) out = poly(Z,degree=3,raw=TRUE);
  if (K==4)  {out = poly(Z,degree=3,raw=TRUE); for (j in 1:ncol(Z)) out = cbind(out,Z[,j]^4);}
  if (K==4.5) out = poly(Z,degree=4,raw=TRUE);
  if (K==5)  {out = poly(Z,degree=4,raw=TRUE); for (j in 1:ncol(Z)) out = cbind(out,Z[,j]^5);}
  if (K==5.5) out = poly(Z,degree=5,raw=TRUE);
  if (K>=6)  {out = poly(Z,degree=5,raw=TRUE); for (k in 6:K) for (j in 1:ncol(Z)) out = cbind(out,Z[,j]^k);}
  ## RETURN POLYNOMIAL BASIS
  return(out)
}

# DGP
data_gen <- function(n) {
  X <- matrix((runif(n*5)-0.5)*2, ncol=5)
  V <- matrix(rnorm(n), ncol=1)
  G <- matrix(exp(diag(X %*% t(X))), ncol=1)
  E <- matrix(0.3637899 * (1 + diag(X %*% t(X))) * V, ncol=1)
  U <- matrix(rnorm(n), ncol=1)
  tt <- matrix(sqrt(diag(X %*% t(X))) + U >= 0, ncol=1) * 1
  Y <- matrix(tt + G + E, ncol=1)
  return(list(Y=Y, X=X, tt=tt))
} 

# cross validation function
K.CV <- function(tt, X, Y) {
  temp <- rep(NaN, nK)
  for (k in 1:nK) {
    X.pol <- cbind(1, tt, gen.P(X, K[k]))
    X.Q   <- qr.Q(qr(X.pol))
    XX <- X.Q %*% t(X.Q)
    Y.hat <- XX %*% Y
    W <- diag(XX)
    temp[k] <- mean(((Y-Y.hat) / (1-W))^2)
  }
  return(which.min(temp))
}

theta.hat2 <- rep(NaN, M)
se.hat2    <- theta.hat2
K.hat2     <- theta.hat2

set.seed(123)
ptm <- proc.time()
for (m in 1:M) {
  data <- data_gen(n)
  X <- data$X; Y <- data$Y; tt <- data$tt
  k.opt <- K.CV(tt, X, Y)
  X.pol <- cbind(1, gen.P(X, K[k.opt]))
  X.Q   <- qr.Q(qr(X.pol))
  MP    <- diag(rep(1,n)) - X.Q %*% t(X.Q)
  Y.M   <- MP %*% Y
  tt.M  <- MP %*% tt
  theta.hat2[m] <- (t(tt.M) %*% Y.M) / (t(tt.M) %*% tt.M)
  Sigma         <- diag((as.numeric((Y.M - tt.M*theta.hat[m, k])))^2)
  se.hat2[m]    <- sqrt(t(tt.M) %*% Sigma %*% tt.M) / (t(tt.M) %*% tt.M)
  K.hat2[m]     <- K.r[k.opt]
}
proc.time() - ptm

# summary of the cross validation 
table(K.hat2)
# estimator
summary(theta.hat2)
sd(theta.hat2)
summary(se.hat2)
sd(se.hat2)

par(mfrow=c(1,2))
hist(theta.hat2, freq=FALSE, xlab="theta-hat", ylab="", main="")
lines(c(mean(theta.hat2), mean(theta.hat2)), c(-1, 20), col="red", lwd=3)
hist(se.hat2, freq=FALSE, xlab="s.e.", ylab="", main="")
lines(c(mean(se.hat2), mean(se.hat2)), c(-1, 80), col="red", lwd=3)

par(mfrow=c(1,2))
CI.l <- theta.hat2 - 1.96 * se.hat2
CI.r <- theta.hat2 + 1.96 * se.hat2
# rejection rate
mean(1 < CI.l | 1 > CI.r)
plot(1:M, CI.l, type="l", ylim=c(0,2), xlab="simulations", ylab="CI")
lines(1:M, CI.r)
abline(1, 0, col="red", lwd=2)

temp <- sort(CI.l, index.return=TRUE)
CI.l <- temp$x
CI.r <- CI.r[temp$ix]
plot(1:M, CI.l, type="l", ylim=c(0,2), xlab="simulations", ylab="CI")
lines(1:M, CI.r)
abline(1, 0, col="red", lwd=2)