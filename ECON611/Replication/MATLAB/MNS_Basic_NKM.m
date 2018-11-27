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

nimpdat 	= 15;	% No. of quarters in impulse response. 

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


% Coefficents for the MA process for Forward Guidance
r1 = 0; 
r2 = 0; 
r3 = 0; 
r4 = 0; 
r5 = 1; 


% ======================================================================= %
% Steady State
% ======================================================================= %


% ======================================================================= %
% Coefficient Matrix
% ======================================================================= %

nlead 	= 1;  	% Number of leads in system 
nlag 	= 1;   	% Number of lags in system 

xnum    = 9;
neq 	= xnum;

% ======================================================================= %
% Position Counters. 
% ======================================================================= %

Xpos		= 1;
Pipos  		= 2;
Rpos        = 3;        % Deviation of real rate from natural rate!
MA1pos   	= 4;		% State variables for MA process (see equation above)
MA2pos   	= 5;	
MA3pos   	= 6;		
MA4pos   	= 7;		
MA5pos   	= 8;
eps_RN_pos  = 9;        % Shock to the real rate

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
Rzero	    = colcount + Rpos;
MA1zero   	= colcount + MA1pos;
MA2zero   	= colcount + MA2pos;
MA3zero   	= colcount + MA3pos;
MA4zero   	= colcount + MA4pos;
MA5zero   	= colcount + MA5pos;
eps_RN_zero	= colcount + eps_RN_pos;

% =========================================================
% Indicators for lead coefficients for each variable: 

colcount  = collead;

Xlead		= colcount + Xpos; 
Pilead		= colcount + Pipos;
Rlead	    = colcount + Rpos;
MA1lead   	= colcount + MA1pos;
MA2lead   	= colcount + MA2pos;
MA3lead   	= colcount + MA3pos;
MA4lead   	= colcount + MA4pos;
MA5lead   	= colcount + MA5pos;
eps_RN_lead	= colcount + eps_RN_pos;

% =========================================================
% Indicators for lag coefficients for each variable: 

colcount  = collag;

Xlag		= colcount + Xpos; 
Pilag		= colcount + Pipos;
Rlag	    = colcount + Rpos;
MA1lag   	= colcount + MA1pos;
MA2lag   	= colcount + MA2pos;
MA3lag   	= colcount + MA3pos;
MA4lag   	= colcount + MA4pos;
MA5lag   	= colcount + MA5pos;
eps_RN_lag	= colcount + eps_RN_pos;

                                    
% ==================================================================================== %
% Coefficient Matrix
% ==================================================================================== %

% Number of coefficients per equation: 
ncoef   = neq*(nlag+nlead+1);

cof     = zeros(neq,ncoef);         % Coef matrix. Each row is an equation 

% IS Curve (in terms of the deviation of real rate from natural rate)

cof(1,Xzero)        = -1;
cof(1,Xlead)        =  1;
cof(1,Rzero)        = -s;


% NKPC 

cof(2,Pizero)       = -1;
cof(2,Pilead)       = beta;
cof(2,Xzero)        = kappa;


% Process for R
cof(3,Rzero)        =  -1;
cof(3,MA1zero)      = r1;
cof(3,MA2zero)      = r2;
cof(3,MA3zero)      = r3;
cof(3,MA4zero)      = r4;
cof(3,MA5zero)      = r5;


% MA Process for RN shock
cof(4,MA1zero)      = -1;
cof(4,eps_RN_zero)  =  1;
cof(5,MA2zero)      = -1;
cof(5,MA1lag)       =  1;
cof(6,MA3zero)      = -1;
cof(6,MA2lag)       =  1;
cof(7,MA4zero)      = -1;
cof(7,MA3lag)       =  1;
cof(8,MA5zero)      = -1;
cof(8,MA4lag)       =  1;


% SHOCK: RN

cof(9,eps_RN_zero)  =  1;



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

shock(eps_RN_pos,1) = -1;

y           = BVAR*shock;                   % initial value 
DATA        = zeros(nimpdat,neq); 
DATA(1,:) 	= y';                           % store initial value 

for t = 2:nimpdat;                          % loop through periods 
 	  y = AVAR*y;
 	  DATA(t,:)=y';
end;

time = 0:(nimpdat-1); 


figure(1) 
plot(time, DATA(:,Rpos),'-r');
hold on; 
plot(time, DATA(:,Xpos), 'b');
%plot(time, DATA(:,Pipos),'-k');
hold off; 
legend('Real interest rate','Output gap');
%legend('Real interest rate','Output gap','Inflation');
xlabel('Horizon in quarters');
ylabel('Percentage points');
ylim([-1.5 1.5])
saveas(figure(1),'IR_1_R','epsc');


%% ====================================================================== %
% ======================================================================= %
%                        Inflation Response
% ======================================================================= %
% ======================================================================= %
%{
% NEED TO THINK ABOUT THIS MORE

inflation = zeros(25,1);
betavec   = zeros(25,1);

for j=0:24
    betavec(j,1)=beta^j;
end

for j=1:25
    inflation(j,1)= sum(betavec(1:j),1);
end
%}


