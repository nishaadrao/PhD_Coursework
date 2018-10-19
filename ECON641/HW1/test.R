## ECON641: PROBLEM SET 1
## Q7: QUANTITATIVE ANALYSIS IN EK MODEL
## Anirudh Yadav 
## 10/18/2018

######################################################################
# Load packages, clear workspace
######################################################################
rm(list = ls())             #clear workspace
library(foreach)            #for looping
library(dplyr)              #for data manipulation
library(data.table)         #for data manipulation
library(ggplot2)            #for pretty plots
library(Matrix)             #fast matrix calcs
options(scipen = 999)       #forces R to use normal numbers instead of scientific notation


######################################################################
# Load WIOD data
######################################################################
data <- read.csv('PhD_Coursework/ECON641/HW1/wiot00_row_apr12.csv',stringsAsFactors=FALSE)

# Remove uneeded rows and columns
data <- data[-c(1,2,3,5,1441:nrow(data)),-c(1,2,4,1440:ncol(data))]