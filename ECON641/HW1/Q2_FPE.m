alpha=0.5;
betaA=0.25;
betaB=0.75;
K_h  = 180;
K_f  = 20;
L_h  = 5;
L_f  = 45;
K    = K_h + K_f;
L    = L_h + L_f;

LintA = (alpha*(1-betaA))/(alpha*(1-betaA)+(1-alpha)*(1-betaB));
LintB = L - LintA;
KintA = (alpha*betaA)/(alpha*betaA+(1-alpha)*betaB);
KintB = K - KintA;