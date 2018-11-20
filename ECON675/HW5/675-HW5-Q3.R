## ECON675: ASSIGNMENT 5
## Q3: WEAK INSTRUMENTS -- EMPIRICAL STUDIES
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
library(boot)               #for bootstrapping
library(mvtnorm)            #for MVN stuff
library(AER)                #for IV regressions
options(scipen = 999)       #forces R to use normal numbers instead of scientific notation

######################################################################
# Input data
###################################################################### 
ak <- fread('PhD_Coursework/ECON675/HW5/Angrist_Krueger.csv')


######################################################################
# [3.1] Angrist-Krueger
######################################################################

# make YOB dummies
ak[, .N, "YoB_ld"]
for(year_i in unique(ak$YoB_ld)){
  
  ak[ ,temp := 0]
  ak[YoB_ld == year_i ,temp := 1]
  setnames(ak, "temp", paste0("d_YOB_ld_", year_i))
  
}
# get a list of all year dummies but one. Exclude the proper one to match coeffs 
year_dummies <- setdiff(grep("d_YOB", colnames(ak), value = TRUE), "d_YOB_ld_0")

# make QoB dummies 
for(qob_i in unique(ak$QoB)){
  
  ak[ ,temp := 0]
  ak[QoB == qob_i ,temp := 1]
  setnames(ak, "temp", paste0("d_QoB_", qob_i))
  
}

# get qob dummy list.  Exclude the proper one to match coeffs 
qob_dummies <- setdiff(grep("d_QoB", colnames(ak), value = TRUE), "d_QoB_1")

# make cross variables of year dummies and qob
#note there is almost certainly a better way to do this but here we are 
inter_list <- NULL
for(d_year in year_dummies){
  
  for(d_qob in qob_dummies){
    
    ak[, temp := get(d_qob)*get(d_year)]
    setnames(ak, "temp", paste0(d_year, "X", d_qob))
    inter_list<- c(inter_list, paste0(d_year, "X", d_qob))
  }
}


# standard controls 
# (i) race, (ii) marrital status, (iii) SMSA, (iv) dummies for
# region, and (iv) dummies for YoB ld.
std_cont <- c("non_white","married", "SMSA", 
              "ENOCENT","ESOCENT", "MIDATL", 
              "MT", "NEWENG", "SOATL", "WNOCENT",
              "WSOCENT",  year_dummies) # get year dummies but leave one out 

# save extra controls 
extra_cont <- c("age_q", "age_sq")

#================#
# ==== ols 1 ====
#================#

# make the formula 
ols1_form <- as.formula(paste0("l_w_wage~educ +", paste(std_cont, collapse = " + ")))

# run ols 
out_ols1 <- data.table(tidy(lm(ols1_form, data = ak)))

# keep what I need 
out_ols1 <- out_ols1[term %chin% c("educ"), c("term", "estimate", "std.error")]
out_ols1[, model := "OLS 1"]

#================#
# ==== OlS 2 ====
#================#

# make the formula 
ols2_form <- as.formula(paste0("l_w_wage~educ +", paste(std_cont, collapse = " + "), " + ", paste0(extra_cont, collapse = " + ")))

#run ols
out_ols2 <- data.table(tidy(lm(ols2_form, data = ak)))

# keep what I need 
out_ols2 <- out_ols2[term %chin% c("educ"), c("term", "estimate", "std.error")]
out_ols2[, model := "OLS 2"]

#===============#
# ==== 2sls ====
#===============#

# write this part as a function so I can use it in 3.2 
# ACTUALLY, im gonna use different faster function but this is fine as a function too
wrap_2sls <- function(in_data){
  
  #=================#
  # ==== 2sls 1 ====
  #=================#
  
  iv_form <- as.formula(paste0("l_w_wage~educ +", paste(std_cont, collapse = " + "), 
                               "| ", 
                               paste(std_cont, collapse = " + "), " + ", paste0(inter_list, collapse = " + ")))
  iv_reg1 <- data.table(tidy(ivreg(iv_form , data = in_data)))
  
  # keep what I need 
  iv_reg1 <- iv_reg1[term %chin% c("educ"), c("term", "estimate", "std.error")]
  iv_reg1[, model := "2sls 1"]
  
  #=================#
  # ==== 2sls 2 ====
  #=================#
  
  
  iv_form2 <- as.formula(paste0("l_w_wage~educ +", paste(std_cont, collapse = " + "), "+", paste0(extra_cont, collapse = " + "),
                                "| ", 
                                paste(std_cont, collapse = " + "), 
                                " + ", paste0(inter_list, collapse = " + "), 
                                "+", paste0(extra_cont, collapse = " + ")))
  iv_reg2 <- data.table(tidy(ivreg(iv_form2 , data = in_data)))
  
  # keep what I need 
  iv_reg2 <- iv_reg2[term %chin% c("educ"), c("term", "estimate", "std.error")]
  iv_reg2[, model := "2sls 2"]
  
  # stack 2sls 
  out_2sls <- rbind(iv_reg1, iv_reg2)
  
  return(out_2sls)
  
}#end 2sls function 

# run function 
ak_2sls <- wrap_2sls(ak)

#========================#
# ==== output tables ====
#========================#

output_3.1 <- rbind(out_ols1, out_ols2, ak_2sls)
setcolorder(output_3.1, c("model", "term", "estimate", "std.error"))

# out put it 

print(xtable(output_3.1, type = "latex"), 
      file = paste0("C://Users/Nmath_000/Documents/Code/courses/econ 675/PS_5_tex/q3.1_table.tex"),
      include.rownames = FALSE,
      floating = FALSE)



#================#
# ==== Q 3.2 ====
#================#


#===================================#
# ==== Fast 2sls function ==========
#===================================#

fast_2sls <- function(in_data){
  
  #===============#
  # ==== reg1 ====
  #===============#
  
  # make x z and y matrices
  y <- as.matrix(in_data[, l_w_wage])
  x <- as.matrix(in_data[, educ])
  cont <- as.matrix(in_data[, c( std_cont, 'const'), with = FALSE])
  z <- as.matrix(in_data[, c(inter_list, std_cont, "const"), with = FALSE])
  
  first_stage_fit <-  z%*%Matrix::solve(Matrix::crossprod(z))%*%(Matrix::crossprod(z, x))
  
  # make x' matrix 
  x_prime <- cbind(first_stage_fit, cont)
  
  
  form_2nd <-  Matrix::solve(Matrix::crossprod(x_prime))%*%(Matrix::crossprod(x_prime, y))
  
  reg1 <- data.table( term = "educ", estimate = form_2nd[1,1], model = "2sls 1")
  
  #===============#
  # ==== reg2 ====
  #===============#
  
  cont <- as.matrix(in_data[, c( std_cont, extra_cont, 'const'), with = FALSE])
  z <- as.matrix(in_data[, c(inter_list, std_cont, extra_cont, "const"), with = FALSE])
  
  first_stage_fit <-  z%*%Matrix::solve(Matrix::crossprod(z))%*%(Matrix::crossprod(z, x))
  
  # make x' matrix 
  x_prime <- cbind(first_stage_fit, cont)
  
  
  form_2nd <-  Matrix::solve(Matrix::crossprod(x_prime))%*%(Matrix::crossprod(x_prime, y))
  
  reg2 <- data.table( term = "educ", estimate = form_2nd[1,1], model = "2sls 2")
  
  # stack results and retur n
  out_results <- rbind(reg1, reg2)
}


#=========================#
# ==== run simulation ====
#=========================#

# copy data for permutation 
ak_perm <- copy(ak)

# add constant 
ak_perm[, const := 1]

# write a function so I can parallel this shiz  
sim_warper <- function(sim_i, in_data = ak_perm ){
  
  # get random sampel 
  perm <- sample(c(1:nrow(in_data)))
  
  # purmute data 
  in_data[, QoB := QoB[perm]]
  
  # clear out dummy variables 
  in_data <- in_data[, -c(grep("d_QoB", colnames(in_data), value = TRUE), inter_list), with = FALSE]
  
  # redo dummy vars 
  for(qob_i in unique(in_data$QoB)){
    
    in_data[ ,temp := 0]
    in_data[QoB == qob_i ,temp := 1]
    setnames(in_data, "temp", paste0("d_QoB_", qob_i))
    
  }
  # recalculate interactions
  inter_list <- NULL
  for(d_year in year_dummies){
    
    for(d_qob in qob_dummies){
      
      in_data[, temp := get(d_qob)*get(d_year)]
      setnames(in_data, "temp", paste0(d_year, "X", d_qob))
      inter_list<- c(inter_list, paste0(d_year, "X", d_qob))
    }
  }
  
  # run 2sls funciton on new data 
  ak_2sls_i <- fast_2sls(in_data)
  
  # add simulation 
  ak_2sls_i[, sim := sim_i]
  
  # retunrn it 
  return(ak_2sls_i)
  
} # end funciton 

# time this sucker 
start_time <- Sys.time()

# parallel setup
cl <- makeCluster(4, type = "PSOCK")
registerDoParallel(cl)

# run simulations in parallel
output_list <- foreach(sim = 1 : 5000,
                       .inorder = FALSE,
                       .packages = "data.table",
                       .options.multicore = list(preschedule = FALSE, cleanup = 9)) %dopar% sim_warper(sim_i = sim)

# stop clusters 
stopCluster(cl)

# check time 
end_time <- Sys.time()

# print time 
print(paste0(round(as.numeric(end_time - start_time, units = "mins"), 3), " minutes to run"))


#==========================#
# ==== organize output ====
#==========================#

# stack data 
sim_res3.2 <- rbindlist(output_list)

# make table 
output3.2 <- sim_res3.2[, list(mean = mean(estimate), std.dev = sd(estimate)), "model"]
