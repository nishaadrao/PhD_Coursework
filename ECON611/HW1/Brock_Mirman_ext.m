% ======================================================================= %
% Brock_Mirman_ext.m
%
% A. Yadav 8/15/2018
% ======================================================================= %

clc; clear all; close all;

nimpdat 	= 200;	% No. of quarters in impulse response. 

% ======================================================================= %
% 			     Parameters 
% ======================================================================= %

% annual rates

r           = .01;		% quarterly discount rate (4 pct annual)
beta        = 1/(1+r); 
deltabar    = .025;     % quarterly depreciation rate (10 pct annual)

% discrete time parameters

s           = .2;       % intertemporal substitution

alpha       = .3;       % capital share in production function

rho_Z       = .95;       % Persistence for income


% Shock Variances

sig_Z       = .01; 
sig_C       = .01; 
sig_D       = .005;

% ======================================================================= %
% Steady State
% ======================================================================= %

K = (alpha/(r+deltabar))^(1/(1-alpha));      
C = K^alpha - deltabar*K;     
Z = 1;      

% ======================================================================= %
% Coefficient Matrix
% ======================================================================= %

nlead 	= 1;  	% Number of leads in system 
nlag 	= 1;   	% Number of lags in system 

xnum    = 6;
neq 	= xnum;

% ======================================================================= %
% Position Counters. 
% ======================================================================= %

Cpos		= 1;
Kpos  		= 2;
Zpos		= 3;
deltapos	= 4;
eps_C_pos   = 5; 
eps_Z_pos   = 6; 

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
Zzero		= colcount + Zpos;
deltazero	= colcount + deltapos;
eps_C_zero	= colcount + eps_C_pos;
eps_Z_zero	= colcount + eps_Z_pos;

% =========================================================
% Indicators for lead coefficients for each variable: 

colcount  = collead;

Clead		= colcount + Cpos; 
Klead		= colcount + Kpos;
Zlead		= colcount + Zpos;
deltalead	= colcount + deltapos;
eps_C_lead	= colcount + eps_C_pos;
eps_Z_lead	= colcount + eps_Z_pos;

% =========================================================
% Indicators for lag coefficients for each variable: 

colcount  = collag;

Clag		= colcount + Cpos; 
Klag		= colcount + Kpos;
Zlag		= colcount + Zpos;
deltalag	= colcount + deltapos;
eps_C_lag	= colcount + eps_C_pos;
eps_Z_lag	= colcount + eps_Z_pos;

                                    
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
cof(1,deltalead)   = beta;

% Resource constraint
% NOTE: here's the key trick: we need to write the LOMK with K_{t+1} as K_t
% and K_t as K_{t-1}. This is similar to what we do in Dynare! 
% Hat tip: TR!
% 

cof(2,Kzero)        = -1;
cof(2,deltazero)    = -1;
cof(2,Klag)        = (1+r);
cof(2,Zzero)        = K^(alpha-1);
cof(2,Czero)        = -C/K;

% Productivity Process

cof(3,Zzero)        = -1;
cof(3,Zlag)         = rho_Z;
cof(3,eps_Z_zero)   = 1;


% SHOCK: Delta

cof(4,deltazero)    = 1;

% SHOCK: C

cof(5, eps_C_zero)  = 1;


% SHOCK: Z

cof(6,eps_Z_zero)   = 1;


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
%                           IMPULSE RESPONSES
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

% MU shock

shock = zeros(neq,1);                       % Shock vector

shock(eps_C_pos,1) = 1;

y           = BVAR*shock;                   % initial value 
DATA        = zeros(nimpdat,neq); 
DATA(1,:) 	= y';                           % store initial value 

for t = 2:nimpdat;                          % loop through periods 
 	  y = AVAR*y;
 	  DATA(t,:)=y';
end;

time = 1:nimpdat; 


figure(2) 
%plot(time, DATA(:,eps_C_pos),'-r');
hold on; 
plot(time, DATA(:,Zpos), 'b');
plot(time, DATA(:,Cpos),'-k');
plot(time, DATA(:,Kpos),'-g');
hold off; 
%legend('e-C','Z','C','K'); 
legend('Z','C','K');
saveas(figure(2),'IR_1_MU','epsc');


% Depreciation shock

shock = zeros(neq,1);                       % Shock vector

shock(deltapos,1) = 1;

y           = BVAR*shock;                   % initial value 
DATA        = zeros(nimpdat,neq); 
DATA(1,:) 	= y';                           % store initial value 

for t = 2:nimpdat;                          % loop through periods 
 	  y = AVAR*y;
 	  DATA(t,:)=y';
end;

time = 1:nimpdat; 


figure(3) 
%plot(time, DATA(:,deltapos),'-r');
hold on; 
plot(time, DATA(:,Zpos), 'b');
plot(time, DATA(:,Cpos),'-k');
plot(time, DATA(:,Kpos),'-g');
hold off; 
%legend('d','Z','C','K'); 
legend('Z','C','K'); 
saveas(figure(3),'IR_1_D','epsc');

%}

%% ====================================================================== %
% ======================================================================= %
%                           SIMULATIONS
% ======================================================================= %
% ======================================================================= %

%% SIMULATE SHOCK PROCESSES 

% These shocks produce `nice' plots of simulated data
deltashock=normrnd(0,sig_D,[1,nimpdat]);
cshock=normrnd(0,sig_C,[1,nimpdat]);
zshock=normrnd(0,sig_Z,[1,nimpdat]);
shock = zeros(neq,nimpdat);
shock([4:6],:)=[deltashock;cshock;zshock];

% Quick plot of shocks
%{
time=1:nimpdat;
plot(time,shock(4,:)',time,shock(5,:)',time,shock(6,:)')
%}

% SIMULATE DATA

y           = zeros(6,1);             % initial value (ie starting at ss)
DATA        = zeros(nimpdat,neq);   

for t = 1:nimpdat                          % loop through periods 
 	  y = AVAR*y+BVAR*shock(:,t);
 	  DATA(t,:)=y';
end

% Quick plot of data

time=1:nimpdat;
plot(time,DATA(:,1),time,DATA(:,2),time,DATA(:,3),time,zeros(1,nimpdat),'-k')
legend('C','K','Z');
%}

%% SIMPLE REGRESSION ANALYSIS

% Get OLS coefficient for the regression of C_t on delta_t

X = DATA(:,4);
Y = DATA(:,1);

betahat = (X'*Y)/(X'*X);

%% ====================================================================== %
% ======================================================================= %
%                       ANALYTICAL VARIANCE OF Y
% ======================================================================= %
% ======================================================================= %
 

% Construct variance matrix for the shocks
omegaE          = zeros(neq,neq);
omegaE(4,4)     = sig_D;
omegaE(5,5)     = sig_C;
omegaE(6,6)     = sig_Z;

% Compute quadratic form with BVAR and vectorize using reshape
BOB             = BVAR*omegaE*BVAR';
BOBvec          = reshape(BOB,[neq*neq,1]);

% Compute analytical var-cov matrix of Y using CH's trick method
Vvec            = inv(eye(neq*neq)- kron(AVAR,AVAR))*BOBvec;
V               = reshape(Vvec,[neq,neq]);

% Compute asymptotic beta
beta0            = V(deltapos,Cpos)/V(deltapos,deltapos);


