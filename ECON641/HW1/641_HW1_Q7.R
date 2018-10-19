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
# Construct bilateral intermeidate trade by country pair
######################################################################
data <- read.csv('PhD_Coursework/ECON641/HW1/wiot00_row_apr12.csv',stringsAsFactors=FALSE)

# Remove uneeded rows and columns
data <- data[-c(1,2,3,5,1441:nrow(data)),-c(1,2,4,1440:ncol(data))]

# Get vector of country names
country.names  <- as.character(data[,1])
country.names  <- country.names[2:length(country.names)]
country.unique <- unique(country.names)

# Remove country names from rows
data  <- data[-1,]

# Remove commas, convert to numbers
data  <- apply(data[,-1], 2, function(y) as.numeric(gsub(",", "", y)))
data  <- as.data.frame(data)

# Add country back as a grouping variable
data  <- cbind(country.names,data)


# Add column names
colnames(data) <- c("country",country.names)

# Convert to data table for easy manipulation
dt             <- as.data.table(data)

# Make vector of industry indicies
ind <- c(paste0("c", c(1:35)))

# Rename cols with industry suffixes 
colnames(dt) <- c(colnames(dt)[1], paste0(colnames(dt)[2:ncol(dt)], "_", ind))


# Sum columns by country
col.sum <- dt[, lapply(.SD, sum), by = "country"]

# Sum rows, by country
# NOTES:
for (i in 1:length(country.unique)){
  col.sum[, temp := rowSums(.SD), .SDcols = grep(paste0(country.unique[i],"_"), colnames(col.sum))]
  setnames(col.sum, "temp", country.unique[i])
}


# subset the data table
intermediate.trade = col.sum[,country.unique,with=FALSE]
intermediate.trade[,supplier:=country.unique]

# reorder the cols
setcolorder(intermediate.trade,c("supplier",country.unique))

# Remove large data
rm(data,dt,col.sum)


######################################################################
# Construct bilateral final goods trade by country pair
######################################################################
data <- read.csv('PhD_Coursework/ECON641/HW1/wiot00_row_apr12.csv',stringsAsFactors=FALSE)

# Remove uneeded rows and columns
data <- data[-c(1,2,3,5,1441:nrow(data)),-c(1,2,4:1439,ncol(data))]

# Get vector of country names
country.names.l  <- as.character(data[,1])
country.names.l  <- country.names.l[2:length(country.names.l)]
country.unique   <- unique(country.names.l)

# Get a short vector of country names
country.names.s  <- as.character(data[1,])
country.names.s  <- country.names.s[2:length(country.names.s)]

# Remove country names row
data  <- data[-1,]

# Remove commas, convert to numbers
data  <- apply(data[,-1], 2, function(y) as.numeric(gsub(",", "", y)))
data  <- as.data.frame(data)

# Add country back as a grouping variable
data  <- cbind(country.names.l,data)


# Add column names
colnames(data) <- c("country",country.names.s)

# Convert to data table for easy manipulation
dt             <- as.data.table(data)

# Make vector of final good indicies
ind <- c(paste0("c", c(1:5)))

# Rename cols with industry suffixes 
colnames(dt) <- c(colnames(dt)[1], paste0(colnames(dt)[2:ncol(dt)], "_", ind))


# Sum columns by country
col.sum <- dt[, lapply(.SD, sum), by = "country"]

# Sum rows, by country
for (i in 1:length(country.unique)){
  col.sum[, temp := rowSums(.SD), .SDcols = grep(paste0(country.unique[i],"_"), colnames(col.sum))]
  setnames(col.sum, "temp", country.unique[i])
}


# subset the data table
final.trade = col.sum[,country.unique,with=FALSE]
final.trade[,supplier:=country.unique]

# reorder the cols
setcolorder(final.trade,c("supplier",country.unique))

# Remove large data
rm(data,dt,col.sum)




