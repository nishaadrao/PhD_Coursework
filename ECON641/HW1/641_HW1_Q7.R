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
data <- read.csv('PhD_Coursework/ECON641/HW1/wiot00_row_apr12_play.csv',stringsAsFactors=FALSE)

# Remove uneeded rows and columns
data <- data[-c(1,2,3,5),-c(1,2,4)]

# Get vector of country names
country.names  <- as.character(data[,1])
country.names  <- country.names[2:length(country.names)]
country.unique <- unique(country.names)

# Remove country names from columns
data  <- data[-1,]

# Remove commas, convert to numbers
data  <- apply(data[,-1], 1, function(y) as.numeric(gsub(",", "", y)))
data  <- as.data.frame(data)

# Add country as a variable
data  <- cbind(country.names,data)


# Add column names
  colnames(data) <- c("country",country.names)
  dt    <- as.data.table(data)
  
  # Make vector of industry indicies
  ind <- c(paste0("c", c(1:35)))

  # Rename cols with industry suffixes 
  colnames(dt) <- c(colnames(dt)[1], paste0(colnames(dt)[2:ncol(dt)], "_", ind))

  
# Sum columns by country
col.sum <- dt[, lapply(.SD, sum), by = "country"]

# Sum rows - NEED TO WRITE A LOOP OVER COUNTRY NAMES
col.sum[, AUS := rowSums(.SD), .SDcols = grep("AUS_", colnames(col.sum))]
col.sum[, AUT := rowSums(.SD), .SDcols = grep("AUT_", colnames(col.sum))]

# subset the data table
intermediate.trade = col.sum[,country.unique,with=FALSE]

