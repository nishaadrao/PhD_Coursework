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
data[,sic1:=cut(sic,c(0,999,1999,3999,4999,5999,6999,8999,9999),labels=c("Ag","MinCon","Man","Tran","WRTrade","Fin","Serv","Pub"))]

  

