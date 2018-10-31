% ======================================================================= %
% MNS_Basic_NKM.m
%
% This file replicates the results from Section I in MNS (2016)
% using a standard "3-equation NKM"
%
% A. Yadav 10/31/2018
% ======================================================================= %
addpath ./AIM_subroutines;
clc; clear all; close all;

nimpdat 	= 40;	% No. of quarters in impulse response. 

% ======================================================================= %
% CALIBRATION
% ======================================================================= %

% annual rates

r           = .01;		% quarterly discount rate (4 pct annual)
beta        = 1/(1+r); 

% discrete time parameters [CHECK THESE]

alpha       =.27;              % If =0 we have a CRS production technology. Else it's decreasing returns to scale (see model equation 5).
theta       =0.698;            % Measure of price stickiness. If =0 then prices are flexible.
lambda      =0.154;            % or alternatively derived endogneously through lambda=(theta^(-1))*(1-theta)*(1-beta*theta)*(1-alpha)/(1-alpha+alpha*epsilon).
phi         =1;                % Elasticity of labor supply.

s           = 1;               % intertemporal substitution

kappa       = lambda*(s+(phi+alpha)/(1-alpha)); % Slope of NKPC


% ======================================================================= %
% Steady State
% ======================================================================= %


% ======================================================================= %
% Coefficient Matrix
% ======================================================================= %

nlead 	= 1;  	% Number of leads in system 
nlag 	= 1;   	% Number of lags in system 

xnum    = 5;
neq 	= xnum;

% ======================================================================= %
% Position Counters. 
% ======================================================================= %

Xpos		= 1;
Pipos  		= 2;
Ipos		= 3;
RNpos       = 4; 
eps_RN_pos  = 5; 

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

Xzero		= colcount + Xpos; 
Pizero		= colcount + Pipos;
Izero		= colcount + Ipos;
RNzero	    = colcount + RNpos;
eps_RN_zero	= colcount + eps_RN_pos;

% =========================================================
% Indicators for lead coefficients for each variable: 

colcount  = collead;

Xlead		= colcount + Xpos; 
Pilead		= colcount + Pipos;
Ilead		= colcount + Ipos;
RNlead	    = colcount + RNpos;
eps_RN_lead	= colcount + eps_RN_pos;

% =========================================================
% Indicators for lag coefficients for each variable: 

colcount  = collag;

Xlag		= colcount + Xpos; 
Pilag		= colcount + Pipos;
Ilag		= colcount + Ipos;
RNlag	    = colcount + RNpos;
eps_RN_lag	= colcount + eps_RN_pos;

                                    
% ==================================================================================== %
% Coefficient Matrix
% ==================================================================================== %

% Number of coefficients per equation: 
ncoef   = neq*(nlag+nlead+1);

cof     = zeros(neq,ncoef);         % Coef matrix. Each row is an equation 

% IS Curve

cof(1,Xzero)        = -1;
cof(1,Xlead)        =  1;
cof(1,Izero)        = -s;
cof(1,Ilead)        =  s;
cof(1,RNzero)       =  s;


% NKPC 

cof(2,Pizero)       = -1;
cof(2,Pilead)       = beta;
cof(2,Xzero)        = kappa;


% Policy rule

cof(3,Izero)        = -1;
cof(3,Pilead)       =  1;
cof(3,RNzero)       =  1;
cof(3,eps_RN_zero)  =  1;


% Process for RN
cof(4,RNzero)       =  1;


% SHOCK: RN

cof(5,eps_RN_zero)  =  1;



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



% Shock to natural rate

shock = zeros(neq,1);                       % Shock vector

shock(eps_RN_pos,1) = 1;

y           = BVAR*shock;                   % initial value 
DATA        = zeros(nimpdat,neq); 
DATA(1,:) 	= y';                           % store initial value 

for t = 2:nimpdat;                          % loop through periods 
 	  y = AVAR*y;
 	  DATA(t,:)=y';
end;

time = 1:nimpdat; 


figure(1) 
plot(time, DATA(:,eps_RN_pos),'-r');
hold on; 
plot(time, DATA(:,Xpos), 'b');
plot(time, DATA(:,Pipos),'-k');
hold off; 
legend('e-R','X','Pi'); 
saveas(figure(1),'IR_1_Z','epsc');

