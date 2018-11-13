## ECON675: ASSIGNMENT 4
## Q2: ESTIMATING AVERAGE TREATMENT EFFECTS
## Anirudh Yadav 
## 11/06/2018

######################################################################
# Load packages, clear workspace
######################################################################
rm(list = ls())             #clear workspace
library(foreach)            #for looping
library(data.table)         #for data manipulation
library(Matrix)             #fast matrix calcs
library(ggplot2)            #for pretty plots
library(sandwich)           #for variance-covariance estimation 
library(xtable)             #for latex tables
library(boot)               #for bootstrapping
library(CausalGAM)          #for computing ATEs
options(scipen = 999)       #forces R to use normal numbers instead of scientific notation

######################################################################
# Input data, add covariates and subset data
######################################################################
data <- as.data.table(read.csv('PhD_Coursework/ECON675/HW4/LaLonde_all.csv'))

data[,log.re74:=log(re74+1)]
data = data[,log.re75:=log(re75+1)]
data = data[,age.sq:=age^2]
data = data[,educ.sq:=educ^2]
data = data[,age.cu:=age^3]
data = data[,black.u74:=black*u74]
data = data[,educ.logre74:=educ*log.re74]

# subset data for LaLonde control only
X.ll = data[treat==1 | treat==0]
Y.ll = data[treat==1 | treat==0,.(re78)]

# subset data for PSID control only
X.ps = data[treat==1 | treat==2]
Y.ps = data[treat==1 | treat==2,.(re78)]

# Recode treatment indicate in PSID control dataset (recode 2's as 0's)
X.ps = X.ps[,treat:=as.numeric(treat==1)]


######################################################################
# Create covariate sets
######################################################################
X.ll.0  = X.ll[,.(treat)]
X.ps.0  = X.ps[,.(treat)]

X.ll.A  = X.ll[,-c("age.sq","educ.sq","age.cu","black.u74","educ.logre74","u74","u75","re78","re74","re75")]
X.ps.A  = X.ps[,-c("age.sq","educ.sq","age.cu","black.u74","educ.logre74","u74","u75","re78","re74","re75")]

X.ll.B  =  X.ll[,-c("age.cu","black.u74","educ.logre74","re78","re74","re75")]
X.ps.B  =  X.ps[,-c("age.cu","black.u74","educ.logre74","re78","re74","re75")]

X.ll.C  =  X.ll[,-c("re78","re74","re75")]
X.ps.C  =  X.ps[,-c("re78","re74","re75")]

######################################################################
# [1] Difference in means
######################################################################
dmeans.ll = lm(as.matrix(Y.ll)~as.matrix(X.ll.0))
dmeans.ps = lm(as.matrix(Y.ps)~as.matrix(X.ps.0))

# Compute robust standard errors
dmeans.ll.se    = sqrt(diag(vcovHC(dmeans.ll, type = "HC1")))
dmeans.ps.se    = sqrt(diag(vcovHC(dmeans.ps, type = "HC1")))

# Compute 95% CIs
dmeans.ll.lower = dmeans.ll$coefficients - 1.96*dmeans.ll.se
dmeans.ll.upper = dmeans.ll$coefficients + 1.96*dmeans.ll.se

dmeans.ps.lower = dmeans.ps$coefficients - 1.96*dmeans.ps.se
dmeans.ps.upper = dmeans.ps$coefficients + 1.96*dmeans.ps.se

# Put results together
dmeans.ll.results = cbind(dmeans.ll$coefficients,dmeans.ll.se,dmeans.ll.lower,dmeans.ll.upper)
dmeans.ps.results = cbind(dmeans.ps$coefficients,dmeans.ps.se,dmeans.ps.lower,dmeans.ps.upper)


######################################################################
# [2] OLS
######################################################################
# Compute OLS coefficients
ols.ll.A = lm(as.matrix(Y.ll)~as.matrix(X.ll.A))
ols.ps.A = lm(as.matrix(Y.ps)~as.matrix(X.ps.A))
ols.ll.B = lm(as.matrix(Y.ll)~as.matrix(X.ll.B))
ols.ps.B = lm(as.matrix(Y.ps)~as.matrix(X.ps.B))
ols.ll.C = lm(as.matrix(Y.ll)~as.matrix(X.ll.C))
ols.ps.C = lm(as.matrix(Y.ps)~as.matrix(X.ps.C))

# Compute robust standard errors
ols.ll.se.A    = sqrt(diag(vcovHC(ols.ll.A, type = "HC1")))
ols.ps.se.A    = sqrt(diag(vcovHC(ols.ps.A, type = "HC1")))
ols.ll.se.B    = sqrt(diag(vcovHC(ols.ll.B, type = "HC1")))
ols.ps.se.B    = sqrt(diag(vcovHC(ols.ps.B, type = "HC1")))
ols.ll.se.C    = sqrt(diag(vcovHC(ols.ll.C, type = "HC1")))
ols.ps.se.C    = sqrt(diag(vcovHC(ols.ps.C, type = "HC1")))

# Compute 95% CIs
ols.ll.lower.A = ols.ll.A$coefficients - 1.96*ols.ll.se.A
ols.ll.upper.A = ols.ll.A$coefficients + 1.96*ols.ll.se.A
ols.ps.lower.A = ols.ps.A$coefficients - 1.96*ols.ps.se.A
ols.ps.upper.A = ols.ps.A$coefficients + 1.96*ols.ps.se.A
ols.ll.lower.B = ols.ll.B$coefficients - 1.96*ols.ll.se.B
ols.ll.upper.B = ols.ll.B$coefficients + 1.96*ols.ll.se.B
ols.ps.lower.B = ols.ps.B$coefficients - 1.96*ols.ps.se.B
ols.ps.upper.B = ols.ps.B$coefficients + 1.96*ols.ps.se.B
ols.ll.lower.C = ols.ll.C$coefficients - 1.96*ols.ll.se.C
ols.ll.upper.C = ols.ll.C$coefficients + 1.96*ols.ll.se.C
ols.ps.lower.C = ols.ps.C$coefficients - 1.96*ols.ps.se.C
ols.ps.upper.C = ols.ps.C$coefficients + 1.96*ols.ps.se.C

# Put treatment effect results together 
ols.ll.results = cbind(c(ols.ll.A$coefficients[2],ols.ll.B$coefficients[2],ols.ll.C$coefficients[2]),c(ols.ll.se.A[2],ols.ll.se.B[2],ols.ll.se.C[2]),c(ols.ll.lower.A[2],ols.ll.lower.B[2],ols.ll.lower.C[2]),c(ols.ll.upper.A[2],ols.ll.upper.B[2],ols.ll.upper.C[2]))
ols.ps.results = cbind(c(ols.ps.A$coefficients[2],ols.ps.B$coefficients[2],ols.ps.C$coefficients[2]),c(ols.ps.se.A[2],ols.ps.se.B[2],ols.ps.se.C[2]),c(ols.ps.lower.A[2],ols.ps.lower.B[2],ols.ps.lower.C[2]),c(ols.ps.upper.A[2],ols.ps.upper.B[2],ols.ps.upper.C[2]))


######################################################################
# [3.A] Regression Imputation, covariate set A
######################################################################

# Subset outcome data for imputation
Y.treat        = data[treat==1,.(re78)]
Y.control.ll   = data[treat==0,.(re78)]
Y.control.ps   = data[treat==2,.(re78)]

# Subset covariates for imputation
X.treat.A       = data[treat==1,-c("age.sq","educ.sq","age.cu","black.u74","educ.logre74","u74","u75","re78","re74","re75","treat")]
X.control.ll.A  = data[treat==0,-c("age.sq","educ.sq","age.cu","black.u74","educ.logre74","u74","u75","re78","re74","re75","treat")]
X.control.ps.A  = data[treat==2,-c("age.sq","educ.sq","age.cu","black.u74","educ.logre74","u74","u75","re78","re74","re75","treat")]

# Get OLS coefficients for imputation
ols.treat.A          = lm(as.matrix(Y.treat)~as.matrix(X.treat.A))
ols.control.ll.A     = lm(as.matrix(Y.control.ll)~as.matrix(X.control.ll.A))
ols.control.ps.A     = lm(as.matrix(Y.control.ps)~as.matrix(X.control.ps.A))

# I need to add constants to the X's to compute imputed treatment effects,
# Then reorder so const is the first variable
X.treat.A[,const:=1]
setcolorder(X.treat.A,c("const"))
X.control.ll.A[,const:=1]
setcolorder(X.control.ll.A,c("const"))
X.control.ps.A[,const:=1]
setcolorder(X.control.ps.A,c("const"))

# Impute `individual treatment effects`
tvec.ri.treat.ll.A      = as.matrix(X.treat.A)%*%(as.vector(ols.treat.A$coefficients)-as.vector(ols.control.ll.A$coefficients))
tvec.ri.treat.ps.A      = as.matrix(X.treat.A)%*%(as.vector(ols.treat.A$coefficients)-as.vector(ols.control.ps.A$coefficients))

tvec.ri.control.ll.A    = as.matrix(X.control.ll.A)%*%(as.vector(ols.treat.A$coefficients)-as.vector(ols.control.ll.A$coefficients))  
tvec.ri.control.ps.A    = as.matrix(X.control.ps.A)%*%(as.vector(ols.treat.A$coefficients)-as.vector(ols.control.ps.A$coefficients))  

# Compute ATEs
ate.ri.ll.A       = mean(c(tvec.ri.treat.ll.A,tvec.ri.control.ll.A))
ate.ri.ps.A       = mean(c(tvec.ri.treat.ps.A,tvec.ri.control.ps.A))

# Compute ATT
att.ri.A          = mean(tvec.ri.treat.ll.A)

######################################################################
# [3.B] Regression Imputation, covariate set B
######################################################################

# Subset covariates for imputation
X.treat.B       = data[treat==1,-c("age.cu","black.u74","educ.logre74","re78","re74","re75","treat")]
X.control.ll.B  = data[treat==0,-c("age.cu","black.u74","educ.logre74","re78","re74","re75","treat")]
X.control.ps.B  = data[treat==2,-c("age.cu","black.u74","educ.logre74","re78","re74","re75","treat")]

# Get OLS coefficients for imputation
ols.treat.B          = lm(as.matrix(Y.treat)~as.matrix(X.treat.B))
ols.control.ll.B     = lm(as.matrix(Y.control.ll)~as.matrix(X.control.ll.B))
ols.control.ps.B     = lm(as.matrix(Y.control.ps)~as.matrix(X.control.ps.B))

# I need to add constants to the X's to compute imputed treatment effects,
# Then reorder so const is the first variable
X.treat.B[,const:=1]
setcolorder(X.treat.B,c("const"))
X.control.ll.B[,const:=1]
setcolorder(X.control.ll.B,c("const"))
X.control.ps.B[,const:=1]
setcolorder(X.control.ps.B,c("const"))

# Impute `individual treatment effects`
tvec.ri.treat.ll.B      = as.matrix(X.treat.B)%*%(as.vector(ols.treat.B$coefficients)-as.vector(ols.control.ll.B$coefficients))
tvec.ri.treat.ps.B      = as.matrix(X.treat.B)%*%(as.vector(ols.treat.B$coefficients)-as.vector(ols.control.ps.B$coefficients))

tvec.ri.control.ll.B    = as.matrix(X.control.ll.B)%*%(as.vector(ols.treat.B$coefficients)-as.vector(ols.control.ll.B$coefficients))  
tvec.ri.control.ps.B    = as.matrix(X.control.ps.B)%*%(as.vector(ols.treat.B$coefficients)-as.vector(ols.control.ps.B$coefficients))  

# Compute ATEs
ate.ri.ll.B       = mean(c(tvec.ri.treat.ll.B,tvec.ri.control.ll.B))
ate.ri.ps.B       = mean(c(tvec.ri.treat.ps.B,tvec.ri.control.ps.B))

# Compute ATT
att.ri.B          = mean(tvec.ri.treat.ll.B)

######################################################################
# [3.C] Regression Imputation, covariate set C
######################################################################

# Subset covariates for imputation
X.treat.C       = data[treat==1,-c("re78","re74","re75","treat")]
X.control.ll.C  = data[treat==0,-c("re78","re74","re75","treat")]
X.control.ps.C  = data[treat==2,-c("re78","re74","re75","treat")]

# Get OLS coefficients for imputation
ols.treat.C          = lm(as.matrix(Y.treat)~as.matrix(X.treat.C))
ols.control.ll.C     = lm(as.matrix(Y.control.ll)~as.matrix(X.control.ll.C))
ols.control.ps.C     = lm(as.matrix(Y.control.ps)~as.matrix(X.control.ps.C))

# I need to add constants to the X's to compute imputed treatment effects,
# Then reorder so const is the first variable
X.treat.C[,const:=1]
setcolorder(X.treat.C,c("const"))
X.control.ll.C[,const:=1]
setcolorder(X.control.ll.C,c("const"))
X.control.ps.C[,const:=1]
setcolorder(X.control.ps.C,c("const"))

# Impute `individual treatment effects`
tvec.ri.treat.ll.C      = as.matrix(X.treat.C)%*%(as.vector(ols.treat.C$coefficients)-as.vector(ols.control.ll.C$coefficients))
tvec.ri.treat.ps.C      = as.matrix(X.treat.C)%*%(as.vector(ols.treat.C$coefficients)-as.vector(ols.control.ps.C$coefficients))

tvec.ri.control.ll.C    = as.matrix(X.control.ll.C)%*%(as.vector(ols.treat.C$coefficients)-as.vector(ols.control.ll.C$coefficients))  
tvec.ri.control.ps.C    = as.matrix(X.control.ps.C)%*%(as.vector(ols.treat.C$coefficients)-as.vector(ols.control.ps.C$coefficients))  

# Compute ATEs
ate.ri.ll.C       = mean(c(tvec.ri.treat.ll.C,tvec.ri.control.ll.C))
ate.ri.ps.C       = mean(c(tvec.ri.treat.ps.C,tvec.ri.control.ps.C))

# Compute ATT
att.ri.C          = mean(tvec.ri.treat.ll.C)

######################################################################
# Compute propensity scores for each sample and model
######################################################################

# Generate treatment outcome variables
T.ll = data[treat==1|treat==0,.(treat)]
T.ps = data[treat==1|treat==2,.(treat)]

#Recode 2's to 0's for PSID sample
T.ps = T.ps[,treat:=as.numeric(treat==1)]

# Get propensity scores using logit regression
prop.ll.A = glm(as.matrix(T.ll) ~ as.matrix(X.ll.A[,-c("treat")]),family = "binomial")
prop.ll.B = glm(as.matrix(T.ll) ~ as.matrix(X.ll.B[,-c("treat")]),family = "binomial")
prop.ll.C = glm(as.matrix(T.ll) ~ as.matrix(X.ll.C[,-c("treat")]),family = "binomial")

prop.ps.A = glm(as.matrix(T.ps) ~ as.matrix(X.ps.A[,-c("treat")]),family = "binomial")
prop.ps.B = glm(as.matrix(T.ps) ~ as.matrix(X.ps.B[,-c("treat")]),family = "binomial")
prop.ps.C = glm(as.matrix(T.ps) ~ as.matrix(X.ps.C[,-c("treat")]),family = "binomial")

# Add prop scores to the data matrices for easy computing of treatment effects
X.ll.ipw = X.ll
X.ll.ipw[,ps.A:=prop.ll.A$fitted.values]
X.ll.ipw[,ps.B:=prop.ll.B$fitted.values]
X.ll.ipw[,ps.C:=prop.ll.C$fitted.values]

X.ps.ipw = X.ps
X.ps.ipw[,ps.A:=prop.ps.A$fitted.values]
X.ps.ipw[,ps.B:=prop.ps.B$fitted.values]
X.ps.ipw[,ps.C:=prop.ps.C$fitted.values]

######################################################################
# [4.A] Inverse Probability Weighting, Lalonde control
######################################################################

# Create variables for computing ATEs 
X.ll.ipw[,t1.A:=treat*re78/ps.A]
X.ll.ipw[,t0.A:=(1-treat)*re78/(1-ps.A)]
X.ll.ipw[,t1.B:=treat*re78/ps.B]
X.ll.ipw[,t0.B:=(1-treat)*re78/(1-ps.B)]
X.ll.ipw[,t1.C:=treat*re78/ps.C]
X.ll.ipw[,t0.C:=(1-treat)*re78/(1-ps.C)]

# Compute proportion of treated respondents
p.ll           = mean(X.ll[,treat])

# Create additional variables for computing ATTs
X.ll.ipw[,t1.att:=treat*re78/p.ll]
X.ll.ipw[,t0.A2:=(1-treat)*re78/(1-ps.A)*(ps.A/p.ll)]
X.ll.ipw[,t0.B2:=(1-treat)*re78/(1-ps.B)*(ps.B/p.ll)]
X.ll.ipw[,t0.C2:=(1-treat)*re78/(1-ps.C)*(ps.C/p.ll)]

# Compute ATEs
ate.ipw.ll.A  = mean(X.ll.ipw[,t1.A])-mean(X.ll.ipw[,t0.A])
ate.ipw.ll.B  = mean(X.ll.ipw[,t1.B])-mean(X.ll.ipw[,t0.B])
ate.ipw.ll.C  = mean(X.ll.ipw[,t1.C])-mean(X.ll.ipw[,t0.C])

# Compute ATTs
att.ipw.ll.A  = mean(X.ll.ipw[,t1.att])-mean(X.ll.ipw[,t0.A2])
att.ipw.ll.B  = mean(X.ll.ipw[,t1.att])-mean(X.ll.ipw[,t0.B2])
att.ipw.ll.C  = mean(X.ll.ipw[,t1.att])-mean(X.ll.ipw[,t0.C2])

######################################################################
# [4.B] Inverse Probability Weighting, PSID control
######################################################################

# Create variables for computing ATEs 
X.ps.ipw[,t1.A:=treat*re78/ps.A]
X.ps.ipw[,t0.A:=(1-treat)*re78/(1-ps.A)]
X.ps.ipw[,t1.B:=treat*re78/ps.B]
X.ps.ipw[,t0.B:=(1-treat)*re78/(1-ps.B)]
X.ps.ipw[,t1.C:=treat*re78/ps.C]
X.ps.ipw[,t0.C:=(1-treat)*re78/(1-ps.C)]

# Compute proportion of treated respondents
p.ps           = mean(X.ps[,treat])

# Create additional variables for computing ATTs
X.ps.ipw[,t1.att:=treat*re78/p.ps]
X.ps.ipw[,t0.A2:=(1-treat)*re78/(1-ps.A)*(ps.A/p.ps)]
X.ps.ipw[,t0.B2:=(1-treat)*re78/(1-ps.B)*(ps.B/p.ps)]
X.ps.ipw[,t0.C2:=(1-treat)*re78/(1-ps.C)*(ps.C/p.ps)]

# Compute ATEs
ate.ipw.ps.A  = mean(X.ps.ipw[,t1.A])-mean(X.ps.ipw[,t0.A])
ate.ipw.ps.B  = mean(X.ps.ipw[,t1.B])-mean(X.ps.ipw[,t0.B])
ate.ipw.ps.C  = mean(X.ps.ipw[,t1.C])-mean(X.ps.ipw[,t0.C])

# Compute ATTs
att.ipw.ps.A  = mean(X.ps.ipw[,t1.att])-mean(X.ps.ipw[,t0.A2])
att.ipw.ps.B  = mean(X.ps.ipw[,t1.att])-mean(X.ps.ipw[,t0.B2])
att.ipw.ps.C  = mean(X.ps.ipw[,t1.att])-mean(X.ps.ipw[,t0.C2])


######################################################################
# [4,5] IPW and Doubly Robust using the "CausalGAM" package
######################################################################
# For some reason my manual IPW estimates did not match STATA's
# Furthermore, I'm not sure how exactly to compute the SEs by hand
# Accordingly, I'm going to use the CausalGAM package below,
# These results match STATA's.

# Covariates A, Lalonde control
ATE.ll.A <- estimate.ATE(pscore.formula = treat ~ age + educ + black + hisp + married + nodegr + log.re74 + log.re75,
                        pscore.family = binomial,
                        outcome.formula.t = re78 ~ age + educ + black + hisp + married + nodegr + log.re74 + log.re75,
                        outcome.formula.c = re78 ~ age + educ + black + hisp + married + nodegr + log.re74 + log.re75,
                        outcome.family = gaussian,
                        treatment.var = "treat",
                        data=as.data.frame(X.ll),
                        divby0.action="t",
                        divby0.tol=0.001,
                        var.gam.plot=FALSE,
                        nboot=0
                        ) 

# Covariates B, Lalonde control
ATE.ll.B <- estimate.ATE(pscore.formula = treat ~ age + educ + black + hisp + married + nodegr + log.re74 + log.re75 + age.sq + educ.sq + u74 + u75,
                         pscore.family = binomial,
                         outcome.formula.t = re78 ~ age + educ + black + hisp + married + nodegr + log.re74 + log.re75 + age.sq + educ.sq + u74 + u75,
                         outcome.formula.c = re78 ~ age + educ + black + hisp + married + nodegr + log.re74 + log.re75 + age.sq + educ.sq + u74 + u75,
                         outcome.family = gaussian,
                         treatment.var = "treat",
                         data=as.data.frame(X.ll),
                         divby0.action="t",
                         divby0.tol=0.001,
                         var.gam.plot=FALSE,
                         nboot=0
) 

# Covariates C, Lalonde control
ATE.ll.C <- estimate.ATE(pscore.formula = treat ~ age + educ + black + hisp + married + nodegr + log.re74 + log.re75 + age.sq + educ.sq + u74 + u75 + age.cu + black.u74 + educ.logre74,
                         pscore.family = binomial,
                         outcome.formula.t = re78 ~ age + educ + black + hisp + married + nodegr + log.re74 + log.re75 + age.sq + educ.sq + u74 + u75 + age.cu + black.u74 + educ.logre74,
                         outcome.formula.c = re78 ~ age + educ + black + hisp + married + nodegr + log.re74 + log.re75 + age.sq + educ.sq + u74 + u75 + age.cu + black.u74 + educ.logre74,
                         outcome.family = gaussian,
                         treatment.var = "treat",
                         data=as.data.frame(X.ll),
                         divby0.action="t",
                         divby0.tol=0.001,
                         var.gam.plot=FALSE,
                         nboot=0
)

# Covariates A, PSID control -- FOR SOME REASON THIS BREAKS?????
# ATE.ps.A <- estimate.ATE(pscore.formula = treat ~ age + educ + black + hisp + married + nodegr + log.re74 + log.re75,
#                          pscore.family = binomial,
#                          outcome.formula.t = re78 ~ age + educ + black + hisp + married + nodegr + log.re74 + log.re75,
#                          outcome.formula.c = re78 ~ age + educ + black + hisp + married + nodegr + log.re74 + log.re75,
#                          outcome.family = gaussian,
#                          treatment.var = "treat",
#                          data=as.data.frame(X.ps),
#                          divby0.action="t",
#                          divby0.tol=0.001,
#                          var.gam.plot=FALSE,
#                          nboot=0,
#                          suppress.warnings = FALSE
# )

# Covariates B, PSID control
ATE.ps.B <- estimate.ATE(pscore.formula = treat ~ age + educ + black + hisp + married + nodegr + log.re74 + log.re75 + age.sq + educ.sq + u74 + u75,
                         pscore.family = binomial,
                         outcome.formula.t = re78 ~ age + educ + black + hisp + married + nodegr + log.re74 + log.re75 + age.sq + educ.sq + u74 + u75,
                         outcome.formula.c = re78 ~ age + educ + black + hisp + married + nodegr + log.re74 + log.re75 + age.sq + educ.sq + u74 + u75,
                         outcome.family = gaussian,
                         treatment.var = "treat",
                         data=as.data.frame(X.ps),
                         divby0.action="t",
                         divby0.tol=0.001,
                         var.gam.plot=FALSE,
                         nboot=0
) 

# Covariates C, PSID control
ATE.ps.C <- estimate.ATE(pscore.formula = treat ~ age + educ + black + hisp + married + nodegr + log.re74 + log.re75 + age.sq + educ.sq + u74 + u75 + age.cu + black.u74 + educ.logre74,
                         pscore.family = binomial,
                         outcome.formula.t = re78 ~ age + educ + black + hisp + married + nodegr + log.re74 + log.re75 + age.sq + educ.sq + u74 + u75 + age.cu + black.u74 + educ.logre74,
                         outcome.formula.c = re78 ~ age + educ + black + hisp + married + nodegr + log.re74 + log.re75 + age.sq + educ.sq + u74 + u75 + age.cu + black.u74 + educ.logre74,
                         outcome.family = gaussian,
                         treatment.var = "treat",
                         data=as.data.frame(X.ps),
                         divby0.action="t",
                         divby0.tol=0.001,
                         var.gam.plot=FALSE,
                         nboot=0
)


######################################################################
#  CONSTRUCT TABLE 1
######################################################################

## LALONDE CONTROL

# Mean Diff + OLS results
a  =    rbind(dmeans.ll.results[2,],ols.ll.results)

# Reg imputation results
b1  =   c(ATE.ll.A$ATE.reg.hat,ATE.ll.A$ATE.reg.asymp.SE,ATE.ll.A$ATE.reg.hat-1.96*ATE.ll.A$ATE.reg.asymp.SE,ATE.ll.A$ATE.reg.hat+1.96*ATE.ll.A$ATE.reg.asymp.SE)
b2  =   c(ATE.ll.B$ATE.reg.hat,ATE.ll.B$ATE.reg.asymp.SE,ATE.ll.B$ATE.reg.hat-1.96*ATE.ll.B$ATE.reg.asymp.SE,ATE.ll.B$ATE.reg.hat+1.96*ATE.ll.B$ATE.reg.asymp.SE)
b3  =   c(ATE.ll.C$ATE.reg.hat,ATE.ll.C$ATE.reg.asymp.SE,ATE.ll.C$ATE.reg.hat-1.96*ATE.ll.C$ATE.reg.asymp.SE,ATE.ll.C$ATE.reg.hat+1.96*ATE.ll.C$ATE.reg.asymp.SE)

# IPW results
c1  =   c(ATE.ll.A$ATE.IPW.hat,ATE.ll.A$ATE.IPW.asymp.SE,ATE.ll.A$ATE.IPW.hat-1.96*ATE.ll.A$ATE.IPW.asymp.SE,ATE.ll.A$ATE.IPW.hat+1.96*ATE.ll.A$ATE.IPW.asymp.SE)
c2  =   c(ATE.ll.B$ATE.IPW.hat,ATE.ll.B$ATE.IPW.asymp.SE,ATE.ll.B$ATE.IPW.hat-1.96*ATE.ll.B$ATE.IPW.asymp.SE,ATE.ll.B$ATE.IPW.hat+1.96*ATE.ll.B$ATE.IPW.asymp.SE)
c3  =   c(ATE.ll.C$ATE.IPW.hat,ATE.ll.C$ATE.IPW.asymp.SE,ATE.ll.C$ATE.IPW.hat-1.96*ATE.ll.C$ATE.IPW.asymp.SE,ATE.ll.C$ATE.IPW.hat+1.96*ATE.ll.C$ATE.IPW.asymp.SE)

# Doubly robust results
d1  =   c(ATE.ll.A$ATE.AIPW.hat,ATE.ll.A$ATE.AIPW.asymp.SE,ATE.ll.A$ATE.AIPW.hat-1.96*ATE.ll.A$ATE.AIPW.asymp.SE,ATE.ll.A$ATE.AIPW.hat+1.96*ATE.ll.A$ATE.AIPW.asymp.SE)
d2  =   c(ATE.ll.B$ATE.AIPW.hat,ATE.ll.B$ATE.AIPW.asymp.SE,ATE.ll.B$ATE.AIPW.hat-1.96*ATE.ll.B$ATE.AIPW.asymp.SE,ATE.ll.B$ATE.AIPW.hat+1.96*ATE.ll.B$ATE.AIPW.asymp.SE)
d3  =   c(ATE.ll.C$ATE.AIPW.hat,ATE.ll.C$ATE.AIPW.asymp.SE,ATE.ll.C$ATE.AIPW.hat-1.96*ATE.ll.C$ATE.AIPW.asymp.SE,ATE.ll.C$ATE.AIPW.hat+1.96*ATE.ll.C$ATE.AIPW.asymp.SE)

## PSID control

# Mean Diff + OLS results
e  =    rbind(dmeans.ps.results[2,],ols.ps.results)

# Reg imputation results
f1  =   c(0,0,0,0)
f2  =   c(ATE.ps.B$ATE.reg.hat,ATE.ps.B$ATE.reg.asymp.SE,ATE.ps.B$ATE.reg.hat-1.96*ATE.ps.B$ATE.reg.asymp.SE,ATE.ps.B$ATE.reg.hat+1.96*ATE.ps.B$ATE.reg.asymp.SE)
f3  =   c(ATE.ps.C$ATE.reg.hat,ATE.ps.C$ATE.reg.asymp.SE,ATE.ps.C$ATE.reg.hat-1.96*ATE.ps.C$ATE.reg.asymp.SE,ATE.ps.C$ATE.reg.hat+1.96*ATE.ps.C$ATE.reg.asymp.SE)

# IPW results
g1  =   c(0,0,0,0)
g2  =   c(ATE.ps.B$ATE.IPW.hat,ATE.ps.B$ATE.IPW.asymp.SE,ATE.ps.B$ATE.IPW.hat-1.96*ATE.ps.B$ATE.IPW.asymp.SE,ATE.ps.B$ATE.IPW.hat+1.96*ATE.ps.B$ATE.IPW.asymp.SE)
g3  =   c(ATE.ps.C$ATE.IPW.hat,ATE.ps.C$ATE.IPW.asymp.SE,ATE.ps.C$ATE.IPW.hat-1.96*ATE.ps.C$ATE.IPW.asymp.SE,ATE.ps.C$ATE.IPW.hat+1.96*ATE.ps.C$ATE.IPW.asymp.SE)

# Doubly robust results
h1  =   c(0,0,0,0)
h2  =   c(ATE.ps.B$ATE.AIPW.hat,ATE.ps.B$ATE.AIPW.asymp.SE,ATE.ps.B$ATE.AIPW.hat-1.96*ATE.ps.B$ATE.AIPW.asymp.SE,ATE.ps.B$ATE.AIPW.hat+1.96*ATE.ps.B$ATE.AIPW.asymp.SE)
h3  =   c(ATE.ps.C$ATE.AIPW.hat,ATE.ps.C$ATE.AIPW.asymp.SE,ATE.ps.C$ATE.AIPW.hat-1.96*ATE.ps.C$ATE.AIPW.asymp.SE,ATE.ps.C$ATE.AIPW.hat+1.96*ATE.ps.C$ATE.AIPW.asymp.SE)


## PUT RESULTS TOGETHER
ll.results  = rbind(a,b1,b2,b3,c1,c2,c3,d1,d2,d3)
ps.results  = rbind(e,f1,f2,f3,g1,g2,g3,h1,h2,h3)
ate.results = round(cbind(ll.results,ps.results),2)

# EXPORT RESULTS AS CSV
setwd("/Users/Anirudh/Desktop/GitHub/PhD_Coursework/ECON675/HW4")
write.table(ate.results, file = "Table1_ATE_resultq.csv",row.names=FALSE,col.names=FALSE,sep=",")

