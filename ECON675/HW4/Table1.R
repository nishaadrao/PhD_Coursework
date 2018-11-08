# The following codes in R generate a tex file appropriate for inclusion in a latex document which displays Table 1 
rm(list=ls())
library(Hmisc)

setwd("/Users/Anirudh/Desktop/GitHub/PhD_Coursework/ECON675/HW4")


tableATEq <- as.matrix(read.csv("Table1_ATE_resultq.csv", header=FALSE))
latex(round(tableATEq, 2),
      file=paste0("teffectsATEq.txt"), append=FALSE,table.env=FALSE,center="none",title="",
      n.cgroup=c(4, 4),
      cgroup=c("Experimental Data", "PSID Control"),
      colheads=c("$\\hat{\\tau}$", "s.e.", "C.I.", "", "$\\hat{\\tau}$", "s.e.", "C.I.", ""),
      n.rgroup=c(1, rep(3, 6)),
      rgroup=c("Mean Diff.", "OLS", "Reg. Impute", "IPW", "D. Robust", "N1 Match", "p Match"),
      rowname=c("", rep(c("a", "b", "c"), 8)))

# Similarly, generate "teffectsATTq.txt"
tableATTq <- as.matrix(read.csv("Table1_ATT_resultq.csv", header=FALSE))
latex(round(tableATTq, 2),
      file=paste0("teffectsATTq.txt"), append=FALSE,table.env=FALSE,center="none",title="",
      n.cgroup=c(4, 4),
      cgroup=c("Experimental Data", "PSID Control"),
      colheads=c("$\\hat{\\tau}$", "s.e.", "C.I.", "", "$\\hat{\\tau}$", "s.e.", "C.I.", ""),
      n.rgroup=c(1, rep(3, 6)),
      rgroup=c("Mean Diff.", "OLS", "Reg. Impute", "IPW", "D. Robust", "N1 Match", "p Match"),
      rowname=c("", rep(c("a", "b", "c"), 8)))
