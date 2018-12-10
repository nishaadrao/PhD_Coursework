## ECON641: PROBLEM SET 2
## Q1: THE FIRM SIZE DISTRIBUTION
## Anirudh Yadav 
## 12/07/2018

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
library(readstata13)        #for reading Stata dta files
options(scipen = 999)       #forces R to use normal numbers instead of scientific notation

######################################################################
# Input data & create new variables
######################################################################
setwd("/Users/Anirudh/Desktop/GitHub/PhD_Coursework/ECON641/HW2")
data = as.data.table(read.dta13('PanelAnnual_compustat1980_2015.dta'))

# Create 1 digit SIC codes
data[,sic1:=cut(sic,c(0,999,1999,3999,4999,5999,6999,8999,9999),labels=c("Agriculture","Mining + Construction","Manufacturing","Transport, etc.","Wholesale + Retail Trade","Finance, etc.","Services","Public Admin"))]

# Create additional variables, by fyear
data[,`:=`(rsales=rank(-sale), remp=rank(-emp), logsale=log(sale),logemp=log(emp)), by=fyear]

# Create industry-level rankings
data[,`:=`(rsales.ind=rank(-sale), remp.ind=rank(-emp)), by=.(fyear,sic1)]

# Keep only 2015 & 1985 data
data=data[fyear==2015|fyear==1985]


######################################################################
# Plot firm size distribution in 2015 -- SALES 
######################################################################

# Plot 1 -- full sample
p2015.full = ggplot(data[fyear==2015],aes(x=logsale,y=log(rsales)))+
                geom_point()+
                #ggtitle("Size Distribution of All Compustat Firms in 2015")+
                #theme(plot.title = element_text(hjust = 0.5))+
                xlab("log(Sales)")+ylab("log(Rank)")
ggsave("2015-sales-full.pdf",p2015.full)

# Plot 2 -- top 500 firms
p2015.500 = ggplot(data[fyear==2015 & rsales<=500],aes(x=logsale,y=log(rsales)))+
                geom_point()+
                xlab("log(Sales)")+ylab("log(Rank)")
ggsave("2015-sales-500.pdf",p2015.500)

# Plot 3 -- top 100 firms
p2015.100 = ggplot(data[fyear==2015 & rsales<=100],aes(x=logsale,y=log(rsales)))+
                geom_point()+
                xlab("log(Sales)")+ylab("log(Rank)")
ggsave("2015-sales-100.pdf",p2015.100)

# Plot 4 -- full sample, by industry
p2015.all.ind = ggplot(data[fyear==2015],aes(x=logsale,y=log(rsales.ind)))+
                 geom_point()+
                 xlab("log(Sales)")+ylab("log(Rank)") + facet_wrap(~sic1)

# Plot 5 -- top 100 firms, by industry
p2015.100.ind = ggplot(data[fyear==2015 & rsales.ind<=100],aes(x=logsale,y=log(rsales.ind)))+
                  geom_point()+
                  xlab("log(Sales)")+ylab("log(Rank)") + facet_wrap(~sic1)

p2015.100.ind.ex = ggplot(data[fyear==2015 & rsales.ind<=100 & sic1!="Agriculture" & sic1!="Public Admin"],aes(x=logsale,y=log(rsales.ind)))+
                    geom_point()+
                    xlab("log(Sales)")+ylab("log(Rank)") + facet_wrap(~sic1)

######################################################################
# Power law coefficients -- SALES
######################################################################

# Different firm size cutoffs
cutoffs = c(length(data[fyear==2015,rsales]),500,100)

# Run Gabaix-Ibra OLS regressions
sales.regs = lapply(1:3,function(i) lm(log(rsales-0.5)~logsale,data[fyear==2015 & rsales<=cutoffs[i] & logsale!=-Inf]))

# Get estimates of PL coefficient
sales.betas = sapply(1:3, function(i) sales.regs[[i]]$coef[2])

# Get OLS SEs
sales.olsSE = sapply(1:3, function(i) summary(sales.regs[[i]])$coefficients[,2][2])

# Compute Gabaix-Ibra SEs
sales.gabaixSE = sapply(1:3, function(i) sqrt(2/cutoffs[i])*abs(sales.betas[i]))

# Combind results
sales.results = rbind(sales.betas,sales.gabaixSE)
colnames(sales.results) = c("Full", "Top 500", "Top 100")


######################################################################
# Power law coefficients -- SALES -- BY INDUSTRY
######################################################################

cutoffs.ind = c(length(data[fyear==2015,rsales]),100,50)
ind.vec = c("Agriculture","Mining + Construction","Manufacturing","Transport, etc.","Wholesale + Retail Trade","Finance, etc.","Services","Public Admin")

# Run OLS regs -- spits out a big list of 3 x 8 x 12 (i.e. 24 regressions)
sales.regs.ind = foreach(j=ind.vec,.combine='cbind')%:%
                  foreach(i=cutoffs.ind, .combine = 'c') %do%{
                   lm(log(rsales.ind-0.5)~logsale,data[fyear==2015 & rsales.ind<=i & logsale!=-Inf & sic1==j])
                }

#lapply(1:3,function(i) lm(log(rsales.ind-0.5)~logsale,data[fyear==2015 & rsales.ind<=cutoffs[i] & logsale!=-Inf] & sic1=ind.vec[j]))


