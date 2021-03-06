## ECON641: PROBLEM SET 1
## Q7: QUANTITATIVE ANALYSIS IN EK MODEL
## Anirudh Yadav 
## 10/18/2018

######################################################################
# Load packages, clear workspace
######################################################################
rm(list = ls())             #clear workspace
library(foreach)            #for looping
library(data.table)         #for data manipulation
library(Matrix)             #fast matrix calcs
library(xtable)
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

######################################################################
# Construct total trade flows by country pair
######################################################################
total.trade = intermediate.trade[,-1]+final.trade[,-1]
total.trade[,supplier:=country.unique]
setcolorder(total.trade,c("supplier",country.unique))
total.trade[,total.sp:=colSums(.SD),.SDcols = country.unique]

######################################################################
# Construct ratio of intermediate imports to total imports, by country
######################################################################

N = length(country.unique)

int.share <- function(i){
  out = intermediate.trade[supplier!=country.unique[i],sum(get(country.unique[i]))]/+
    (total.trade[supplier!=country.unique[i],sum(get(country.unique[i]))])
}

intermediate.share = sapply(1:N,function(i) int.share(i))

mean.int.share = mean(intermediate.share)

intermediate.share = as.data.table(round(intermediate.share,2))
intermediate.share[,country:=country.unique]
colnames(intermediate.share)=c("intermediate.share",colnames(intermediate.share)[2])
setcolorder(intermediate.share,c("country","intermediate.share"))


######################################################################
# Construct ratio of trade deficit to total expenditure, by country
######################################################################

total.exp.f <- function(i){
  out = total.trade[supplier==country.unique[i],rowSums(.SD),.SDcols=country.unique[-i]]
}

total.imp.f <- function(i){
  out = total.trade[supplier!=country.unique[i],sum(get(country.unique[i]))]
}

total.exp = sapply(1:N, function(i) total.exp.f(i))
total.imp = sapply(1:N, function(i) total.imp.f(i))
total.sp  = total.trade[,total.sp]

trade.deficits = as.data.table(cbind(total.exp,total.imp,total.sp))
trade.deficits[,deficit:=total.imp-total.exp]
trade.deficits[,ratio:=round(deficit/total.sp,3)]
trade.deficits[,country:=country.unique]
setcolorder(trade.deficits,c("country","total.imp","total.exp","total.sp","ratio"))

mean.def.ratio = trade.deficits[,mean(ratio)]

######################################################################
# Construct bilateral trade shares
######################################################################

# Make matrix of total spending by country
m1=t(matrix(rep(total.trade[,total.sp],2),41,41))

# Get total bilateral trade values as a matrix
m2=as.matrix(total.trade[,country.unique,with=FALSE])

# Get bilateral trade shares
bilateral.shares = as.data.table(m2/m1)

# Round shares to 3 decimal places
bilateral.shares.r = round(bilateral.shares,3)

bilateral.shares.r[,supplier:=country.unique]
setcolorder(bilateral.shares.r,c("supplier",country.unique))

# Remove intermediate matricies
rm(m2)


######################################################################
# 7.2 Counterfactuals
######################################################################
theta       = 8.28
T.hat       = 1
d.hat       = 0.9
d.hat.mat   = matrix(d.hat,41,41)+diag(rep(1-d.hat,41))
nu          = 0.2

# Initial guess of no change
w.hat0      = rep(1,N)

# Get vector of data discrepencies (probably due to value added)
D   =     rowSums(total.trade[,-c(1,43)])-total.trade[,total.sp]

# Enforce initial world GDP normalization
X              = total.trade[,total.sp]/sum(total.trade[,total.sp])

# Write function to update guess of w.hat given an initial guess
update.guess   = function(w.hat0,T.hat=1,d.hat=1){
    
    # Create matrix of initial guesses
    w.hat0.mat  = matrix(w.hat0,41,41)
    
    # Compute counterfactual trade shares
    pi.numerator   = bilateral.shares*(w.hat0.mat*d.hat.mat)^(-theta)
    pi.denominator = t(matrix(colSums(bilateral.shares*(w.hat0.mat*d.hat.mat)^(-theta)),41,41))
    pi.dash        = pi.numerator/pi.denominator
    
    # COMPUTE EXCESS DEMAND VECTOR
    Z                 = rowSums(pi.dash*t(w.hat0.mat)*m1) - w.hat0*total.trade[,total.sp] - D
    
    # UPDATE GUESS USING EXCESS DEMAND
    w.hat.update      = w.hat0*(1+nu*Z/total.trade[,total.sp])
    
    # Renormalize the new guess
    T.w.x             = w.hat.update*X/sum(w.hat.update*X)
    
    # Compute updated guess
    w.hat1            = T.w.x/X

    return(w.hat1)
}   

# LOOOOOOOOOOOOOP until convergence
tol        = 10^(-4) 
maxiter    = 100
converged  = FALSE 
iter       = 0

while (converged == FALSE & iter < maxiter){
  
  w.hat1 = update.guess(w.hat0)
  
  diff   = w.hat1 - w.hat0
  
  converged = all(diff<tol)
  
  w.hat0 = w.hat1
  
  iter   = iter + 1
}

# COMPUTE WELFARE CHANGES
w.hat.mat      = matrix(w.hat1,41,41)
pi.numerator   = bilateral.shares*(w.hat0.mat*d.hat.mat)^(-theta)
pi.denominator = t(matrix(colSums(bilateral.shares*(w.hat0.mat*d.hat.mat)^(-theta)),41,41))
pi.dash        = pi.numerator/pi.denominator

welfare        = (diag(as.matrix(bilateral.shares))/diag(as.matrix(pi.dash)))^(1/theta)

# PUT RESULTS IN A TABLE FOR LATEX
cf.results = as.data.table(cbind(w.hat1,welfare))
cf.results[,Country:=country.unique]
setcolorder(cf.results,c("Country","w.hat1","welfare"))
