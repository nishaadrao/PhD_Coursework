% ======================================================================= %
% Brock_Mirman_ext.m
%
% A. Yadav 8/15/2018
% ======================================================================= %
addpath ./AIM_subroutines;
clc; clear all; close all;

nimpdat 	= 200;	% No. of quarters in impulse response. 

% ======================================================================= %
% 			     Parameters 
% ======================================================================= %

% annual rates

r           = .02;		% quarterly discount rate (4 pct annual)
beta        = 1/(1+r); 
deltabar    = .025;     % quarterly depreciation rate (10 pct annual)

% discrete time parameters

s           = 1;       % intertemporal substitution

alpha       = .35;       % capital share in production function

rho_Z       = .95;       % Persistence for income


% Shock Variances

sig_Z       = .01; 
sig_C       = .002; 
sig_D       = .002;
sig_G       = .002;

% ======================================================================= %
% Steady State
% ======================================================================= %

K = (alpha/(r+deltabar))^(1/(1-alpha));
Y = K^alpha;
I = deltabar*K;
C = K^alpha - deltabar*K;     
Z = 1;      

% ======================================================================= %
% Coefficient Matrix
% ======================================================================= %

nlead 	= 1;  	% Number of leads in system 
nlag 	= 1;   	% Number of lags in system 

xnum    = 8;
neq 	= xnum;

% ======================================================================= %
% Position Counters. 
% ======================================================================= %

Cpos		= 1;
Kpos  		= 2;
Ipos        = 3;
Zpos		= 4;
eps_Z_pos   = 5; 
eps_D_pos	= 6;
eps_C_pos   = 7; 
eps_G_pos   = 8; 

% ==================================================================================== %

colzero = 0+nlag*xnum;      % Position counter for start of contemp. coefs 
                            % i.e. if there is one lag and 7 equations then the "first"
                            % contemporaneous variable is at position 8 (8 = 1*7 + ypos)
collead = 0+nlag*xnum+xnum; % Position counter for start of lead coefs 
                            % i.e. if there is one lag and 7 equations then the "first"
                            % 7 lags are accounted for in nlag*xnum then the contemporaneous
                            % variables take up positions 8 - 14 so that the first lead variable 
                            % is at position 15 (nlag*xnum + xnum + cpos)
collag  = 0;                % Position counter for start of lag coefs  
                            % i.e. these variables start immediately. (place 1)


% =========================================================
% Indicators for contemporanous coefs for each variable: 

colcount  = colzero;

Czero		= colcount + Cpos; 
Kzero		= colcount + Kpos;
Izero       = colcount + Ipos;
Zzero		= colcount + Zpos;
eps_Z_zero	= colcount + eps_Z_pos;
eps_D_zero	= colcount + eps_D_pos;
eps_C_zero	= colcount + eps_C_pos;
eps_G_zero  = colcount + eps_G_pos;

% =========================================================
% Indicators for lead coefficients for each variable: 

colcount  = collead;

Clead		= colcount + Cpos; 
Klead		= colcount + Kpos;
Ilead       = colcount + Ipos;
Zlead		= colcount + Zpos;
eps_Z_lead	= colcount + eps_Z_pos;
eps_D_lead	= colcount + eps_D_pos;
eps_C_lead	= colcount + eps_C_pos;
eps_G_lead  = colcount + eps_G_pos;

% =========================================================
% Indicators for lag coefficients for each variable: 

colcount  = collag;

Clag		= colcount + Cpos; 
Klag		= colcount + Kpos;
Ilag        = colcount + Ipos;
Zlag		= colcount + Zpos;
eps_Z_lag	= colcount + eps_Z_pos;
eps_D_lag	= colcount + eps_D_pos;
eps_C_lag	= colcount + eps_C_pos;
eps_G_lag   = colcount + eps_G_pos;

                                    
% ==================================================================================== %
% Coefficient Matrix
% ==================================================================================== %

% Number of coefficients per equation: 
ncoef   = neq*(nlag+nlead+1);

cof     = zeros(neq,ncoef);         % Coef matrix. Each row is an equation 

% Euler Equation

cof(1,Czero)        = 1/s;
cof(1,eps_C_zero)   = -1;
cof(1,Clead)        = -1/s;
cof(1,eps_C_lead)   = 1;
cof(1,Zlead)        = beta*(r+deltabar);
cof(1,Kzero)        = -beta*(1-alpha)*(r+deltabar);

% Capital accumulation
% NOTE: we need to write the LOMK with K_{t+1} as K_t

cof(2,Kzero)        = -1;
cof(2,Klag)         = (1-deltabar);
cof(2,Izero)        = deltabar;
cof(2,eps_D_zero)     = -1;

% Resource constraint

cof(3,Zzero)        = -1;
cof(3,Klag)         = -alpha;
cof(3,eps_G_zero)   = -1;
cof(3,Izero)        = I/Y;
cof(3,Czero)        = C/Y;

% Productivity Process

cof(4,Zzero)        = -1;
cof(4,Zlag)         = rho_Z;
cof(4,eps_Z_zero)   = 1;

% SHOCK: Z

cof(5,eps_Z_zero)   = 1;

% SHOCK: Delta

cof(6,eps_D_zero)   = 1;

% SHOCK: C

cof(7, eps_C_zero)  = 1;

% SHOCK: Z

cof(8, eps_G_zero)  = 1;


% ==================================================================================== %
%			End coefficient matrix setup
% ==================================================================================== %

% ==================================================================================== %
%		        Begin Solution Algorithm 
% ==================================================================================== %
%                                                                 
%  Solve a linear perfect foresight model using the gauss eig     
%  function to find the invariant subspace associated with the big 
%  roots.  This procedure will fail if the companion matrix is     
%  defective and does not have a linearly independent set of       
%  eigenvectors associated with the big roots.                     
%                                                                  
%  Input arguments:                                                
%                                                                  
%    h         Structural coef matrix (neq,neq*(nlag+1+nlag)).     
%    neq       Number of equations.                                
%    nlag      Number of lags.                                     
%    nlag      Number of lags.                                      
%    condn     lag tolerance used as a condition number test       
%              by numeric_shift and reduced_form.                  
%    uprbnd    Inclusive upper bound for the modulus of roots      
%              allowed in the reduced form.                        
%                                                                  
%  Output arguments:                                               
%                                                                  
%    cofb      Reduced form coefficient matrix (neq,neq*nlag).     
%    scof      Observable Structure                                
%    amat      Companion form matrix                               
%    b         Contemporaneous coefficient matrix                  
%    Model satisfies:                                              
%                                                                  
%    z(t) = amat*z(t-1) + b*e(t)                                   
%                                                                  
%    where the first neq elements of z(t) are the contemporaneous  
%    values of the variables in the model                          
%    and e(t) is the shock vector of conformable dimension         
%                                                                  
%    rts       Roots returned by eig.                              
%    ia        Dimension of companion matrix (number of non-trivial
%              elements in rts).                                   
%    nexact    Number of exact shiftrights.                        
%    nnumeric  Number of numeric shiftrights.                      
%    lgroots   Number of roots greater in modulus than uprbnd.     
%    mcode     Return code: see function aimerr.                   
% ==================================================================================== %


% Use AIM procedure to solve model: 
uprbnd = 1+1e-8;    % Tolerance values for AIM program 
condn = 1e-8;

% ==================================================================================== %
% Run AIM
% ==================================================================================== %

[cofb,rts,ia,nex,nnum,lgrts,mcode] = ...
       aim_eig(cof,neq,nlag,nlead,condn,uprbnd);

scof = obstruct(cof,cofb,neq,nlag,nlead);

% ==================================================================================== %
% need to calculate amat and b
% ==================================================================================== %

s0 = scof(:,(neq*nlag+1):neq*(nlag+1)); 	% Contemp. coefs from obs. structure
amat=zeros(neq*nlag,neq*nlag);   		% Initialize A matrix 
bmat=cofb(1:neq,((nlag-1)*neq+1):nlag*neq);  	% Lag 1 coefficients 
i=2;
while i<=nlag;
  bmat=[bmat cofb(1:neq,((nlag-i)*neq+1):(nlag-i+1)*neq)]; % Lag i coefs 
  i=i+1;
end;
amat(1:neq,:)=bmat;  % Coefs for equations 
if nlag>1;
 amat((length(cofb(:,1))+1):length(amat(:,1)),1:neq*(nlag-1))=eye(neq*(nlag-1));
end;
b = zeros(length(amat(:,1)),length(s0(1,:)));
b(1:length(s0(:,1)),1:length(s0(1,:))) = inv(s0);  % Store coefs 
%b=b(:,shockvec);

AVAR = amat; 
BVAR = b; 

%% ====================================================================== %
% ======================================================================= %
%                      IMPULSE RESPONSES (CHECK)
% ======================================================================= %
% ======================================================================= %
%{
% Productivity shock

shock = zeros(neq,1);                       % Shock vector

shock(eps_Z_pos,1) = 1;

y           = BVAR*shock;                   % initial value 
DATA        = zeros(nimpdat,neq); 
DATA(1,:) 	= y';                           % store initial value 

for t = 2:nimpdat;                          % loop through periods 
 	  y = AVAR*y;
 	  DATA(t,:)=y';
end;

time = 1:nimpdat; 


figure(1) 
plot(time, DATA(:,eps_Z_pos),'-r');
hold on; 
plot(time, DATA(:,Zpos), 'b');
plot(time, DATA(:,Cpos),'-k');
plot(time, DATA(:,Kpos),'-g');
hold off; 
legend('e-Z','Z','C','K'); 
saveas(figure(1),'IR_1_Z','epsc');
%}

%% ====================================================================== %
% ======================================================================= %
%                           RECOVER SHOCKS
% ======================================================================= %
% ======================================================================= %
load('Brock_Mirman_DATA.mat');

% Construct square A and B matrices
A = AVAR(1:4,1:4);
B = BVAR(1:4,5:8);


% Recover shocks from simulated data
shock_rec       = zeros(size(DATA,1),4);   

for t = 2:size(DATA,1)                         
 	  u = inv(B)*(DATA(t,:)'-A*DATA(t-1,:)');
 	  shock_rec(t,:)=u';
end

% Plot shocks
time=1:size(DATA,1);
figure(1)
plot(time,shock_rec(:,1),time,shock_rec(:,2),time,shock_rec(:,3),time,shock_rec(:,4))
legend('e_Z','e_D','e_C','e_G');


%% ====================================================================== %
% ======================================================================= %
%                       SIMULATIONS (CHECK)
% ======================================================================= %
% ======================================================================= %
%% SIMULATE SHOCK PROCESSES 
%{
% These shocks produce `nice' plots of simulated data
e_Z=normrnd(0,sig_Z,[1,nimpdat]);
e_D=normrnd(0,sig_D,[1,nimpdat]);
e_C=normrnd(0,sig_C,[1,nimpdat]);
e_G=normrnd(0,sig_G,[1,nimpdat]);

shock = zeros(neq,nimpdat);
shock([5:8],:)=[e_Z;e_D;e_C;e_G];

% Quick plot of shocks

figure(2)
time=1:nimpdat;
plot(time,shock(5,:)',time,shock(6,:)',time,shock(7,:)',time,shock(8,:)')
legend('e_Z','e_D','e_C','e_G');

% SIMULATE DATA

y           = zeros(neq,1);             % initial value (ie starting at ss)
DATA        = zeros(nimpdat,neq);   

for t = 1:nimpdat                          % loop through periods 
 	  y = AVAR*y+BVAR*shock(:,t);
 	  DATA(t,:)=y';
end

% Quick plot of data

figure(3)
time=1:nimpdat;
plot(time,DATA(:,1),time,DATA(:,2),time,DATA(:,3),time,DATA(:,4),time,zeros(1,nimpdat),'-k')
legend('C','K','I','Z');


% Recover shocks from simulated data
shock_rec_sim       = zeros(nimpdat,neq);   

for t = 2:nimpdat                          
 	  u = inv(BVAR)*(DATA(t,:)'-AVAR*DATA(t-1,:)');
 	  shock_rec_sim(t,:)=u';
end

figure(4)
plot(time,shock_rec_sim(:,5),time,shock_rec_sim(:,6),time,shock_rec_sim(:,7),time,shock_rec_sim(:,8))
legend('e_Z','e_D','e_C','e_G');
%}


%% ====================================================================== %
% ======================================================================= %
%                           LOG LIKELIHOOD
% ======================================================================= %
% ======================================================================= %

% Create var-cov matrix
sigma = diag([sig_Z,sig_D,sig_C,sig_G]);


% Compute MVN pdf values using handwritten function
pdf_vals       = zeros(size(DATA,1),1);

for t = 2:size(DATA,1)                         
 	  p = mymvn(shock_rec(t,:),sigma);
 	  pdf_vals(t,1)=p;
end

% Cross-check using MATLAB's mvnpdf function
pdf_vals_check       = zeros(size(DATA,1),1);

for t = 2:size(DATA,1)                         
 	  q = mvnpdf(shock_rec(t,:),[0,0,0,0],sigma);
 	  pdf_vals_check(t,1)=q;
end


% Compute log-likelihood

log_pdf = log(pdf_vals);

L       = sum(log_pdf(2:size(pdf_vals,1),1));


%% ====================================================================== %
% ======================================================================= %
%                           EXTRA QUESTION 
% ======================================================================= %
% ======================================================================= %

alpha_vec = 0.11:0.01:0.40;
delta_vec = 0.01:0.01:0.10;


[ALPHA,DELTA] = meshgrid(alpha_vec,delta_vec);

LMAT          = zeros(length(delta_vec),length(alpha_vec));

% Run nested loop to compute likelihoods for different combos of alpha and
% delta

for d=1:length(delta_vec)

    deltabar = delta_vec(d);

    for a=1:length(alpha_vec)

        alpha = alpha_vec(a);
        
        LMAT(d,a)=mybrockmirman(alpha,deltabar);

    end

end

%% MAKE SURFACE PLOT
figure(5)
surf(ALPHA,DELTA,LMAT)
title('Log Likelihoodvalues for Different Values of (\alpha,\delta)')
xlabel('\alpha')
ylabel('\delta')
zlabel('L')

% Get alpha/delta combo that maximizes L
[Lmax,Idx] = max(LMAT(:));
[LmaxRow,LmaxCol] = ind2sub(size(LMAT), Idx);
MLE_vec = [alpha_vec(LmaxCol),delta_vec(LmaxRow)];
