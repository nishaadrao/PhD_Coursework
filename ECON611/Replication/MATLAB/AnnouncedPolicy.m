% ======================================================================= %
% BasicRANKchris.m
% For ECON 611, Replication Project
% John S. Olson
% This code performs forward guidance in a super basic model
% Mostly programmed by our Lord and Savior Chris House, God bless Him.
% ======================================================================= %

clc; clear all; close all;

nimpdat = 10;	% No. of quarters in impulse response. 



% ======================================================================= %
% 			     Parameters 
% ======================================================================= %

% Parameters
alpha   	= 0.36;			% share of capital in ouput
beta    	= 0.99;			% discount factor
delta  		= 0.025;		% depreciation of capital
sigmaC		= 1;			% risk aversion consumption
kappa       = .4;           % coeffcient on Y/C in a basic Phillips Curve
rho_r		= 0.7;			% Monetary Policy Smoothing Parameter

% Shock Variances
sig_R = 1;

% Coefficents for the MA process for Forward Guidance
% This process is MA = r1*MA1 + r2*MA2 + r3*MA3 + r4*MA4 + r5*MA5
% For forward guidance we only want the interest rate to drop in period 5
% Hence we only have a non-zero coefficent in period 5.
% See class notes for a tax example where the shock is an immediate tax
% change but agents know it will expire in a few years.
r1 = 0; 
r2 = 0; 
r3 = 0; 
r4 = 0; 
r5 = 1; 



% ======================================================================= %
% Coefficient Matrix
% ======================================================================= %

nlead = 1;  	% Number of leads in system 
nlag = 1;   	% Number of lags in system 

xnum = 9; 
neq = xnum;

% ======================================================================= %
% Position Counters. 
% ======================================================================= %

var_x	= 1;		% Consumption or Income
var_pi	= 2;		% Inflation
var_r	= 3;		% Real Interest Rate
var_MA1	= 4;		% State variables for MA process (see equation above)
var_MA2	= 5;		
var_MA3	= 6;		
var_MA4	= 7;		
var_MA5	= 8;		
e_rPOS  = 9;        % Interest Rate Shock
% ======================================================================= %

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

var_xZero = colcount + var_x;	
var_piZero = colcount + var_pi;
var_rZero = colcount + var_r;	
var_MA1zero = colcount + var_MA1;
var_MA2zero = colcount + var_MA2;
var_MA3zero = colcount + var_MA3;
var_MA4zero = colcount + var_MA4;
var_MA5zero = colcount + var_MA5;
e_rzero = colcount + e_rPOS;     


% =========================================================
% Indicators for lead coefficients for each variable: 

colcount  = collead;

var_xLead = colcount + var_x;	
var_piLead = colcount + var_pi;
var_rLead = colcount + var_r;	
var_MA1lead = colcount + var_MA1;
var_MA2lead = colcount + var_MA2;
var_MA3lead = colcount + var_MA3;
var_MA4lead = colcount + var_MA4;
var_MA5lead = colcount + var_MA5;
e_rlead = colcount + e_rPOS;     


% =========================================================
% Indicators for lag coefficients for each variable: 

colcount  = collag;

var_xLag = colcount + var_x;	
var_piLag = colcount + var_pi;
var_rLag = colcount + var_r;	
var_MA1lag = colcount + var_MA1;
var_MA2lag = colcount + var_MA2;
var_MA3lag = colcount + var_MA3;
var_MA4lag = colcount + var_MA4;
var_MA5lag = colcount + var_MA5;
e_rLag = colcount + e_rPOS;     
                   


% ======================================================================= %
% Coefficient Matrix
% ======================================================================= %

ncoef   = neq*(nlag+nlead+1); % Number of coefficients per equation: 
cof     = zeros(neq,ncoef); % Coef matrix. Each row is an equation 

% Euler Equation
cof(1,var_xLead) = 1;
cof(1,var_xZero) = -1;
cof(1,var_rZero) = -sigmaC;

% Price dynamics 
cof(2,var_piZero) = -1;
cof(2,var_piLead) = beta;
cof(2,var_xZero) = kappa;

% Monetary Policy 
% Note that we are including the state variables for the MA process here
% Because here Fed policy is declaring the states for forward guidiance
% i.e. "This is how interest rates are gonna be for the next 5 periods."
cof(3,var_rZero) = -1;
cof(3,var_MA1zero) = r1;
cof(3,var_MA2zero) = r2;
cof(3,var_MA3zero) = r3;
cof(3,var_MA4zero) = r4;
cof(3,var_MA5zero) = r5;


% Evolution of MA terms 
% In the first period, the announcement of lower rates comes as a shock
% Hence the shock enters in the first period. The remaining terms are 
% passing this information along. Since the agents know the monetary policy
% process, they know the shock won't affect things until the fifth period. 
% But since they know the process they can plan for it.
cof(4,var_MA1zero) = -1;
cof(4,e_rzero) = 1;
cof(5,var_MA2zero) = -1;
cof(5,var_MA1lag) =  1;
cof(6,var_MA3zero) = -1;
cof(6,var_MA2lag) =  1;
cof(7,var_MA4zero) = -1;
cof(7,var_MA3lag) =  1;
cof(8,var_MA5zero) = -1;
cof(8,var_MA4lag) =  1;

% SHOCK: Monetary Policy
cof(9,e_rzero) = 1;

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
% ======================================================================= %

% Use AIM procedure to solve model: 
uprbnd = 1+1e-8;    % Tolerance values for AIM program 
condn = 1e-8;

% ======================================================================= %
% Run AIM
% ======================================================================= %

[cofb,rts,ia,nex,nnum,lgrts,mcode] = ...
       aim_eig(cof,neq,nlag,nlead,condn,uprbnd);

scof = obstruct(cof,cofb,neq,nlag,nlead);

% ======================================================================= %
% need to calculate amat and b
% ======================================================================= %

s0 = scof(:,(neq*nlag+1):neq*(nlag+1)); 	% Contemp. coefs from obs. structure
amat=zeros(neq*nlag,neq*nlag);   		% Initialize A matrix 
bmat=cofb(1:neq,((nlag-1)*neq+1):nlag*neq);  	% Lag 1 coefficients 
i=2;
while i<=nlag
  bmat=[bmat cofb(1:neq,((nlag-i)*neq+1):(nlag-i+1)*neq)]; % Lag i coefs 
  i=i+1;
end
amat(1:neq,:)=bmat;  % Coefs for equations 
if nlag>1
 amat((length(cofb(:,1))+1):length(amat(:,1)),1:neq*(nlag-1))=eye(neq*(nlag-1));
end
b = zeros(length(amat(:,1)),length(s0(1,:)));
b(1:length(s0(:,1)),1:length(s0(1,:))) = inv(s0);  % Store coefs 
%b=b(:,shockvec);

AVAR = amat; 
BVAR = b; 



% ======================================================================= %
%                           IMPULSE RESPONSE
% ======================================================================= %

% Forward Guidance Shock
shock = zeros(neq,1);                       % Shock vector
shock(e_rPOS,1) = 1;

y = BVAR*shock;                   % initial value 
DATA1 = zeros(nimpdat,neq); 
DATA1(1,:) = y';                           % store initial value 

for t = 2:nimpdat                          % loop through periods 
 	  y = AVAR*y;
 	  DATA1(t,:)=y';
end

time = 1:nimpdat; 

figure(1) 
plot(time, DATA1(:,var_x),'LineWidth',3);
hold on; 
plot(time, DATA1(:,var_pi),'LineWidth',3);
plot(time, DATA1(:,var_r),'LineWidth',3);
hold off; 
legend('x', 'pi', 'r'); 


