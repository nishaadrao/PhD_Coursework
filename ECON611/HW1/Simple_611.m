% ======================================================================= %
% Simple_611.m
%
% C. House 9/07/14
% ======================================================================= %

clc; clear all; close all;

nimpdat 	= 20;	% No. of quarters in impulse response. 

% ======================================================================= %
% 			     Parameters 
% ======================================================================= %

% annual rates

r           = .01;		% quarterly discount rate (4 pct annual)
beta        = 1/(1+r); 

% discrete time parameters

s           = .2;       % intertemporal substitution.

rho_Y       = .5;       % Persistence for income 
rho_r       = .5;       % Persistence for real int. rate 

% Shock Variances

sig_Y       = .01; 
sig_R       = .02; 

% ======================================================================= %
% Steady State
% ======================================================================= %

Y = 1;      
C = 1;      
S = 0;      

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
Ypos  		= 2;
Spos		= 3;
rpos		= 4;
eps_Y_pos   = 5; 
eps_r_pos   = 6; 

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
Yzero		= colcount + Ypos;
Szero		= colcount + Spos;
rzero		= colcount + rpos;
eps_Y_zero	= colcount + eps_Y_pos;
eps_r_zero	= colcount + eps_r_pos;

% =========================================================
% Indicators for lead coefficients for each variable: 

colcount  = collead;

Clead		= colcount + Cpos; 
Ylead		= colcount + Ypos;
Slead		= colcount + Spos;
rlead		= colcount + rpos;
eps_Y_lead	= colcount + eps_Y_pos;
eps_r_lead	= colcount + eps_r_pos;

% =========================================================
% Indicators for lag coefficients for each variable: 

colcount  = collag;

Clag		= colcount + Cpos; 
Ylag		= colcount + Ypos;
Slag		= colcount + Spos;
rlag		= colcount + rpos;
eps_Y_lag	= colcount + eps_Y_pos;
eps_r_lag	= colcount + eps_r_pos;

                                    
% ==================================================================================== %
% Coefficient Matrix
% ==================================================================================== %

% Number of coefficients per equation: 
ncoef   = neq*(nlag+nlead+1);

cof     = zeros(neq,ncoef);         % Coef matrix. Each row is an equation 

% Euler Equation

cof(1,Czero)       = -1;
cof(1,Clead)       = 1;
cof(1,rzero)       = -s*beta;


% Flow Budget Constraint

cof(2,Czero)        = -1;
cof(2,Szero)        = -1;
cof(2,Yzero)        = 1;
cof(2,Slag)         = 1+r;


% Income Process

cof(3,Yzero)        = -1;
cof(3,Ylag)         = rho_Y;
cof(3,eps_Y_zero)   = 1;


% Interest Rate Process

cof(4,rzero)        = -1;
cof(4,rlag)         = rho_r;
cof(4,eps_r_zero)   = 1;


% SHOCK: Y

cof(5, eps_Y_zero)  = 1;


% SHOCK: r

cof(6,eps_r_zero)   = 1;


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

% ======================================================================= %
% ======================================================================= %
%                           IMPULSE RESPONSES
% ======================================================================= %
% ======================================================================= %

% Income shock

shock = zeros(neq,1);                       % Shock vector

shock(eps_Y_pos,1) = 1;

y           = BVAR*shock;                   % initial value 
DATA        = zeros(nimpdat,neq); 
DATA(1,:) 	= y';                           % store initial value 

for t = 2:nimpdat;                          % loop through periods 
 	  y = AVAR*y;
 	  DATA(t,:)=y';
end;

time = 1:nimpdat; 


figure(1) 
plot(time, DATA(:,Ypos),'-r');
hold on; 
plot(time, DATA(:,Cpos), 'b');
plot(time, DATA(:,rpos),'-k');
plot(time, DATA(:,Spos),'-m');
plot(time, DATA(:,eps_Y_pos),'-g');
hold off; 
legend('Y', 'C', 'r', 'S','e'); 

%{

% Interest Rate Shock 

shock = zeros(neq,1);                       % Shock vector

shock(eps_r_pos,1) = 1;

y           = BVAR*shock;                   % initial value 
DATA        = zeros(nimpdat,neq); 
DATA(1,:) 	= y';                           % store initial value 

for t = 2:nimpdat;                          % loop through periods 
 	  y = AVAR*y;
 	  DATA(t,:)=y';
end;

time = 1:nimpdat; 


figure(2) 
plot(time, DATA(:,Ypos),'-r');
hold on; 
plot(time, DATA(:,Cpos), 'b');
plot(time, DATA(:,rpos),'-k');
plot(time, DATA(:,Spos),'-m');
hold off; 
legend('Y', 'C', 'r', 'S'); 

%}


