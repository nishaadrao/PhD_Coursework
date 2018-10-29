## ECON641: PROBLEM SET 1
## Q2: HO MODEL
## Anirudh Yadav 
## 10/28/2018


######################################################################
# Calibration
######################################################################
alpha    =0.5
betaA    =0.25
betaB    =0.75
K_h      =180
K_f      =20
L_h      = 5
L_f      = 45
K        = K_h + K_f
L        = L_h + L_f


######################################################################
# [2] FPE set
######################################################################
LintA = ((alpha*(1-betaA))/(alpha*(1-betaA)+(1-alpha)*(1-betaB)))*L
LintB = L - LintA
KintA = ((alpha*betaA)/(alpha*betaA+(1-alpha)*betaB))*K
KintB = K - KintA

x = c(LintA,LintB,LintA+LintB,L_h)
y = c(KintA,KintB,KintA+KintB,K_h)

plot(c(LintA,LintB,LintA+LintB,L_h),c(KintA,KintB,KintA+KintB,K_h),xlim=c(0,50),main="FPE Set",ylim=c(0,200),xlab="L",ylab="K",col=ifelse(x==L_h, "red", "black"))
arrows(0,0,LintA,KintA)
arrows(0,0,LintB,KintB)
arrows(LintB,KintB,LintA+LintB,KintA+KintB)
arrows(LintA,KintA,LintA+LintB,KintA+KintB)

