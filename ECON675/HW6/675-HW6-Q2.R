## ECON675: ASSIGNMENT 5
## Q2: WEAK INSTRUMENTS SIMULATIONS
## Anirudh Yadav 
## 11/19/2018

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
library(rdrobust)           #for RD plots and other stuff
library(rddensity)          #for RD density continuity tests
library(rdlocrand)          #for RD randomization inference
options(scipen = 999)       #forces R to use normal numbers instead of scientific notation

######################################################################
# Input data
######################################################################
setwd("/Users/Anirudh/Desktop/GitHub/PhD_Coursework/ECON675/HW6")
data <- as.data.table(read.csv('HeadStart.csv'))

######################################################################
# [2.1]A RD Plots of pre-intervention mortality
######################################################################

# Evenly-spaced bins, IMSE optimal
rdplot(data[,mort_related_pre],data[,povrate60],p=1,binselect = "es",x.label="povrate60",y.label="mort_related_pre",title="")
dev.copy(pdf,'q2-1-es.pdf')
dev.off()

# Evenly-spaced bins, mimicking variance 
rdplot(data[,mort_related_pre],data[,povrate60],p=1,binselect = "esmv",x.label="povrate60",y.label="mort_related_pre",title="")
dev.copy(pdf,'q2-1-esmv.pdf')
dev.off()

# Quantile-spaced bins, IMSE optimal
rdplot(data[,mort_related_pre],data[,povrate60],p=1,binselect = "qs",x.label="povrate60",y.label="mort_related_pre",title="")
dev.copy(pdf,'q2-1-qs.pdf')
dev.off()

# Quantile-spaced bins, mimicking variance
rdplot(data[,mort_related_pre],data[,povrate60],p=1,binselect = "qsmv",x.label="povrate60",y.label="mort_related_pre",title="")
dev.copy(pdf,'q2-1-qsmv.pdf')
dev.off()

######################################################################
# [2.1]B Formal falsification tests
######################################################################

## Exact binomial tests for different windows around the cutoff

# Vector of windows
h.vec = seq(0.3,1.3,0.2)

# Get running variable
x     = data[,povrate60]

# Number of observations just above and below the cutoff
N.l   = sapply(1:length(h.vec),function(i) sum(x >= -h.vec[i] & x <=0))
N.u   = sapply(1:length(h.vec),function(i) sum(x >= 0 & x <= h.vec[i]))

# Total number of observations in the window
N.t   = N.l + N.u

# Conduct exact binomial tests (p=0.5), where success is treatment and store p-vals
binom.pvals = sapply(1:length(h.vec),function(i) binom.test(N.u[i],N.t[i])$p.value)

# Put results together for latex
binom.results = cbind(h.vec,N.l,N.u,binom.pvals)
xtable(binom.results,digits = c(0,1,0,0,3))

## Continuity in density tests (defaults are triangular kernel, jackknife SEs)
rdtest = rddensity(x)

######################################################################
# [2.2]A Global polynomial regression - constant treatment effect
######################################################################

# Create treatment dummy for regressions
treat=ifelse(data[,povrate60]>=0,1,0)

# Get outcome variable
Y = data[,mort_related_post]

# Generate covariates for polynomial regressions
X.pol = cbind(x,x^2,x^3,x^4,x^5,x^6)

# Run polynomial regressions
global.regs = lapply(0:3,function(i) lm(Y ~ treat + X.pol[,c(1:(3+i))]))

# Get point estimates
global.betas = sapply(1:4,function(i) global.regs[[i]]$coefficients[2])

# Get robust SEs
global.SEs   = sapply(1:4,function(i) sqrt(diag(vcovHC(global.regs[[i]],"HC2")))[2])

# Put results together
global.results = rbind(global.betas, global.SEs)
colnames(global.results) = c(3,4,5,6)
xtable(global.results,digits=c(0,2,2,2,2))

# Plot fitted values and data
temp.rd = rdplot(Y,x,hide=TRUE)
plot(temp.rd$vars_bins$rdplot_mean_x,temp.rd$vars_bins$rdplot_mean_y,pch=20,xlab="povrate60",ylab="mort_related_post")
points(x,global.regs[[2]]$fitted.values,pch=6,col="blue")
abline(v=0)
dev.copy(pdf,'q2-2-const.pdf')
dev.off()

######################################################################
# [2.2]B Global polynomial regression - fully interacted model
######################################################################

# Run fully-interacted polynomial regressions
global.regs.full = lapply(0:3,function(i) lm(Y ~ treat + X.pol[,c(1:(3+i))] + treat*X.pol[,c(1:(3+i))]))

# Get point estimates of treatment effect
global.betas.full = sapply(1:4,function(i) global.regs.full[[i]]$coefficients[2])

# Get robust SEs
global.SEs.full   = sapply(1:4,function(i) sqrt(diag(vcovHC(global.regs.full[[i]],"HC2")))[2])

# Put results together
global.results.full = rbind(global.betas.full, global.SEs.full)
colnames(global.results.full) = c(3,4,5,6)
xtable(global.results.full,digits=c(0,2,2,2,2))

# Plot fitted values + data
temp.rd = rdplot(Y,x,hide=TRUE)
plot(temp.rd$vars_bins$rdplot_mean_x,temp.rd$vars_bins$rdplot_mean_y,pch=20,xlab="povrate60",ylab="mort_related_post")
points(x,global.regs.full[[2]]$fitted.values,pch=6,col="blue")
abline(v=0)
dev.copy(pdf,'q2-2-full.pdf')
dev.off()

######################################################################
# [2.3] Robust local polynomial methods
######################################################################

# MSE-optimal RD treatment effect estimates
rd.regs = lapply(0:2, function(i) rdrobust(Y,x,p=i,all=TRUE))

# Combine results for different polynomial orders
rd.p0   = cbind(rd.regs[[1]]$coef,rd.regs[[1]]$se,rd.regs[[1]]$ci)
rd.p1   = cbind(rd.regs[[2]]$coef,rd.regs[[2]]$se,rd.regs[[2]]$ci)
rd.p2   = cbind(rd.regs[[3]]$coef,rd.regs[[3]]$se,rd.regs[[3]]$ci)

xtable(rd.p0,digits=c(0,2,2,2,2))
xtable(rd.p1,digits=c(0,2,2,2,2))
xtable(rd.p2,digits=c(0,2,2,2,2))

# Placebo tests
rd.rob1 = rdrobust(data[,mort_related_pre],x,p=1,all=TRUE)
rd.rob2 = rdrobust(data[,mort_injury_post],x,p=1,all=TRUE)

rd.rob1.res   = cbind(rd.rob1$coef,rd.rob1$se,rd.rob1$ci)
rd.rob2.res   = cbind(rd.rob2$coef,rd.rob2$se,rd.rob2$ci)

xtable(rd.rob1.res,digits=c(0,2,2,2,2))
xtable(rd.rob2.res,digits=c(0,2,2,2,2))

# Different kernel/bandwiths

######################################################################
# [2.3] Donut holes
######################################################################

# Order data according to running variable
data.ordered = data[order(povrate60)]
x.ordered    = data.ordered[,povrate60]
y.ordered    = data.ordered[,mort_related_post]

# Get indexes of closest counties to the cutoff
i.below        = max(which(x.ordered<=0))
i.above        = min(which(x.ordered>=0))

# Run rdrobust for each donut hole!
rd.donut = lapply(0:9,function(i) rdrobust(y.ordered[-seq(i.below-i,i.above+i,1)],x.ordered[-seq(i.below-i,i.above+i,1)]))

# Get point estimates
rd.donut.results = numeric()
for (t in 1:10){
  ans = rd.donut[[t]]$coef[3]
  rd.donut.results = append(rd.donut.results,ans)
}

######################################################################
# [2.3] Placebo cutoffs
######################################################################

cutoffs = seq(-10,10,2)

rd.cutoffs = lapply(1:length(cutoffs), function(i) rdrobust(Y,x,c=cutoffs[i]))

rd.cutoffs.results = numeric()
for (t in 1:length(cutoffs)){
  ans = rbind(rd.cutoffs[[t]]$coef[3],rd.cutoffs[[t]]$pv[3])
  rd.cutoffs.results = cbind(rd.cutoffs.results,ans)
}


######################################################################
# [2.4] Local randomization inference -- FISHER
######################################################################

# Use defaults to compute recommended window for local randomization
rdwindow = rdwinselect(x,c(data[,mort_related_pre],data[,mort_injury_post]))

# Conduct randomization inference using recommended window
rd.rand.res = rdrandinf(Y,x,wl=rdwindow$window[1],wr=rdwindow$window[2])

######################################################################
# [2.4] Local randomization inference -- NEYMAN
######################################################################

windows = seq(0.8,2.6,0.2)

N.0   = sapply(1:length(windows),function(i) sum(x >= -windows[i] & x <=0))
N.1   = sapply(1:length(windows),function(i) sum(x >= 0 & x <= windows[i]))

T.diffmeans = sapply(1:length(windows), function(i) data[povrate60 >= 0 & povrate60 <=windows[i], mean(mort_related_post)]-data[povrate60 >= -windows[i] & povrate60 <=0, mean(mort_related_post)])

SD.0  = sapply(1:length(windows), function(i) data[povrate60 >= -windows[i] & povrate60 <=0, sd(mort_related_post)] )
SD.1  = sapply(1:length(windows), function(i) data[povrate60 >= 0 & povrate60 <=windows[i], sd(mort_related_post)] )

