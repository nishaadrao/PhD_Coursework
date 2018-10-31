function [ resid , J_L, J, J_P] = ZLB_equ( X, eps )

global Params;


%       1    2    3    4      5     6     7           8       9         10       11
nms = {'Y', 'C', 'N', 'ppi', 'ii', 'S', 'dividend', 'pstar', 'pbarA', 'pbarB', 'wage'};


nnms = length(nms);

for nmi = 1:nnms
    eval([nms{nmi} '_L  = X(nmi);']);
end
nmi_0 = nnms;
for nmi = 1:nnms
    eval([nms{nmi} '  = X(nmi_0+nmi);']);
end
nmi_0 = 2*nnms;
for nmi = 1:nnms
    eval([nms{nmi} '_P  = X(nmi_0+nmi);']);
end    



resid = [...
Y*S - N;   %1
Y - C;     %2
-ii + max( Params.Rbar-1 + Params.phiP*(ppi-1) + eps(1) ,0);       %3
-ppi + ((1-Params.theta)/(1-Params.theta*pstar.^(1/(1-Params.mu)))).^(1-Params.mu);       %4
-pstar + pbarA/pbarB;     %5
-S + (1-Params.theta)*S_L*ppi.^(-Params.mu/(1-Params.mu)) + Params.theta*pstar.^(Params.mu/(1-Params.mu));        %6
-pbarB + Y + Params.beta *  eps(2) *(1-Params.theta)*ppi_P.^(-Params.mu/(1-Params.mu)-1)*pbarB_P;      %7
-pbarA + wage*Params.mu*Y  + Params.beta *  eps(2) *(1-Params.theta)*ppi_P.^(-Params.mu/(1-Params.mu))*pbarA_P;      %8
-dividend +  Y - wage * N;     %9
-Params.psi1*N.^(Params.psi2) + C.^(-Params.sigma)*wage * eps(3);      %10
-C.^(-Params.sigma) + eps(4) *  Params.beta * eps(2) *(1+ii)/ppi_P * C_P.^(-Params.sigma)];     %11

if nargout == 1, return; end


% linearize

iY = 1;
iC = 2;
iN = 3;
ippi = 4;
iii = 5;
iS = 6;
idividend = 7;
ipstar = 8;
ipbarA = 9;
ipbarB = 10;
iwage = 11;

%derivatives wrt lead variables
J_P = zeros(length(resid),nnms);

J_P(7,ippi) = (-Params.mu/(1-Params.mu)-1) * Params.beta *  eps(2) *(1-Params.theta)*ppi_P.^(-Params.mu/(1-Params.mu)-2)*pbarB_P;
J_P(7,ipbarB) = Params.beta *  eps(2) *(1-Params.theta)*ppi_P.^(-Params.mu/(1-Params.mu)-1);

J_P(8,ipbarA) = Params.beta *  eps(2) *(1-Params.theta)*ppi_P.^(-Params.mu/(1-Params.mu));
J_P(8,ippi) = (-Params.mu/(1-Params.mu))*Params.beta *  eps(2) *(1-Params.theta)*ppi_P.^(-Params.mu/(1-Params.mu)-1)*pbarA_P;

J_P(11,ippi) = -eps(4) *Params.beta * eps(2) *(1+ii)/ppi_P^2 * C_P.^(-Params.sigma);
J_P(11,iC) =  (-Params.sigma)*eps(4) *Params.beta * eps(2) *(1+ii)/ppi_P * C_P.^(-Params.sigma-1);


%derivatives wrt lagged variables
J_L = zeros(length(resid),nnms);

J_L(6,iS) = (1-Params.theta)*ppi.^(-Params.mu/(1-Params.mu));





%derivatives wrt current variables
J = zeros(length(resid),nnms);
J(1,iY) = S;
J(1,iS) = Y;
J(1,iN) = -1;

J(2,iY) = 1;
J(2,iC) = -1;

J(3,iii) = -1;
if Params.Rbar-1 + Params.phiP*(ppi-1) + eps(1)  > 0
    J(3,ippi) = Params.phiP;
else
    J(3,ippi) = 0;
end

J(4,ippi) = -1;
J(4,ipstar) = (1-Params.mu) * ((1-Params.theta)/(1-Params.theta*pstar.^(1/(1-Params.mu)))).^(-Params.mu) ...
*(1/(1-Params.mu)) * Params.theta*((1-Params.theta)/(1-Params.theta*pstar.^(1/(1-Params.mu)))^2)*pstar.^(1/(1-Params.mu)-1);

J(5,ipstar) = -1;
J(5,ipbarA) = 1/pbarB;
J(5,ipbarB) = -pbarA/pbarB^2;

J(6,iS) = -1;
J(6,ippi) = (-Params.mu/(1-Params.mu))*(1-Params.theta)*S_L*ppi.^(-Params.mu/(1-Params.mu)-1);
J(6,ipstar) = (Params.mu/(1-Params.mu)) * Params.theta*pstar.^(Params.mu/(1-Params.mu)-1);

J(7,ipbarB) = -1;
J(7,iY) = 1;

J(8,ipbarA) = -1;
J(8,iY) = wage*Params.mu;
J(8,iwage) = Y*Params.mu;

J(9,idividend) = -1;
J(9,iY) = 1;
J(9,iwage) = -N;
J(9,iN) = -wage;

J(10,iN) = -Params.psi1 * Params.psi2 * N.^(Params.psi2-1);
J(10,iC) = (-Params.sigma)*eps(3) * C.^(-Params.sigma-1)*wage;
J(10,iwage) =  eps(3)* C.^(-Params.sigma);

J(11,iC) = -(-Params.sigma)*C.^(-Params.sigma-1);
J(11,iii) = eps(4) *Params.beta * eps(2) /ppi_P * C_P.^(-Params.sigma);


end

