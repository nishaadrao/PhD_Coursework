% ======================================================================= %
% HW1_611.m v2 with "K"~K+1 (a control)
% Simple extension to the Brock Mirman (1972) model
% Tereza Ranosova
% ======================================================================= %
% Q: what is the meaning of coefficients in BVAR to thinks that are never
% shocked?
% Q: EVEN MORE what is the meaning of coefficients picking up a shock to
% "control" (that should never happen?)
% Q: in euler equation, is there a difference if I ignore delta and epsC,
% since those are 0 in expectation?
% For what its worth, I think these two versions are equivalent

clc; clear all; close all;

nimpdat 	= 50;	% No. of quarters in impulse response. 

% ======================================================================= %
% 			     Parameters 
% ======================================================================= %


% discrete time parameters - taken mostly from 607 lectures
r           = .01;	    % quarterly discount rate (4 pct annual)
beta        = 1/(1+r); 
sigma       = .3;       % intertemporal substitution
rhoZ        = .95;       % Persistence for productivity shocks 
delta_bar   = .02;      % ss depreciation
alpha       = .4;       % from RK/Y share?

% Shock Variances - is there any guidence to pick these?(in percentages for
% Z and eps, but absolute for delta)
sig_Z       = .01;     % shock to productivity
sig_D       = .005;     % shock to depreciation (1 a perc point)
sig_Pr      = .01;     % shock to preferences (MUc)

% ======================================================================= %
% Steady State
% ======================================================================= %

K = ((1/beta + delta_bar -1)/alpha)^(1/(alpha-1));      
C = K^alpha - delta_bar*K;      
Z = 1;


% ======================================================================= %
% Coefficient Matrix
% ======================================================================= %

nlead 	= 1;  	% Number of leads in system 
nlag 	= 1;   	% Number of lags in system 

xnum    = 6;    % number of variables. (structural inov as variables!)
neq 	= xnum; % it better be? Well I think it should be - shock_num

% ======================================================================= %
% Position Counters. 
% ======================================================================= %

Cpos		= 1;
Kpos  		= 2; % should not respond to shocks at this moment
Zpos		= 3;
deltapos    = 4;
eps_Z_pos   = 5; 
eps_C_pos   = 6; 

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
eps_Z_zero	= colcount + eps_Z_pos;
eps_C_zero	= colcount + eps_C_pos;

% =========================================================
% Indicators for lead coefficients for each variable: 

colcount  = collead;

Clead		= colcount + Cpos; 
Klead		= colcount + Kpos;
Zlead		= colcount + Zpos;
deltalead	= colcount + deltapos;
eps_Z_lead	= colcount + eps_Z_pos;
eps_C_lead	= colcount + eps_C_pos;

% =========================================================
% Indicators for lag coefficients for each variable: 

colcount  = collag;

Clag		= colcount + Cpos; 
Klag		= colcount + Kpos;
Zlag		= colcount + Zpos;
deltalag	= colcount + deltapos;
eps_Z_lag	= colcount + eps_Z_pos;
eps_C_lag	= colcount + eps_C_pos;
                                
% ==================================================================================== %
% Coefficient Matrix
% ==================================================================================== %

% Number of coefficients per equation: 
ncoef   = neq*(nlag+nlead+1);

cof     = zeros(neq,ncoef);         % Coef matrix. Each row is an equation 

% Euler Equation

cof(1,Czero)       = 1/sigma;
cof(1,Clead)       = -1/sigma;
cof(1,eps_C_zero)  = -1;
%cof(1,eps_C_lead)  = 1;
cof(1,Zlead)       = (1/beta +delta_bar-1)*beta;
cof(1,Kzero)       = (1/beta +delta_bar-1)*beta*(alpha-1);
%cof(1,deltalead)   = -beta;


% Law of motion of kapital - should I write it as Kt=...? (yes, if in lead,
% looks like if it was an expected value, not a determined state
% shouldnt I then as in dynare write product with K_1?

cof(2,Klag)        = 1/beta;
cof(2,Kzero)        = -1;
cof(2,deltazero)    = -1; % this is questionable, but if Kzero=K+1, should be affected by deltazero
cof(2,Zzero)        = K^(alpha-1);
cof(2,Czero)        = -C/K;

% Productivity Process

cof(3,Zzero)        = -1;
cof(3,Zlag)         = rhoZ;
cof(3,eps_Z_zero)   = 1;

%%% order the shocks as they are in the list of variables? Is it important?
% SHOCK: delta
cof(4, deltazero)    = 1;

% SHOCK: Z
cof(5,eps_Z_zero)   = 1;

% SHOCK: Mu(C)
cof(6,eps_C_zero)   = 1;


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
%    z(t) = amat*z(t-1) + b*e(t)     % so in this specitication these are actually "expectation errors I think                              
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

scof = obstruct(cof,cofb,neq,nlag,nlead); % there  is some reshuffling of the variables. the coeffs are put up for some reason

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

AVAR = amat; % propagation (something is wrong in here)
BVAR = b; % initial contemporenous shock

% ======================================================================= %
% ======================================================================= %
%                           IMPULSE RESPONSES
% ======================================================================= %
% ======================================================================= %

%% MU shock

shock = zeros(neq,1);                       % Shock vector
shock(eps_C_pos,1) = 1*sig_Pr;

y           = BVAR*shock;                   % initial value (ie starting at ss)
DATA        = zeros(nimpdat,neq); 
DATA(1,:) 	= y';                           % store initial value 

for t = 2:nimpdat;                          % loop through periods 
 	  y = AVAR*y;
 	  DATA(t,:)=y';
end;

time = 1:nimpdat; 

figure(1) 
plot(time, DATA(:,Cpos),'-r');              % the shocked variable
hold on; 
plot(time, DATA(:,Kpos), 'b');
plot(time, DATA(:,Zpos),'-k');              %even though this should not move?
hold off; 
legend('C', 'K', 'Z'); 


%% Depreciation Shock 

shock = zeros(neq,1);                       % Shock vector

shock(deltapos,1) = 1*sig_D;

y           = BVAR*shock;                   % initial value (ie starting at ss)
DATA        = zeros(nimpdat,neq); 
DATA(1,:) 	= y';                           % store initial value 

for t = 2:nimpdat;                          % loop through periods 
 	  y = AVAR*y;
 	  DATA(t,:)=y';
end;

time = 1:nimpdat; 
figure(2) 
plot(time, DATA(:,Cpos),'-r');              
hold on; 
plot(time, DATA(:,Kpos), 'b');
plot(time, DATA(:,Zpos),'-k');              %even though this should not move?
hold off; 
legend('C', 'K', 'Z');

%% Productiviy Shock 

shock = zeros(neq,1);                       % Shock vector

shock(eps_Z_pos,1) = 1*sig_Z;

y           = BVAR*shock;                   % initial value (ie starting at ss)
DATA        = zeros(nimpdat,neq); 
DATA(1,:) 	= y';                           % store initial value 

for t = 2:nimpdat;                          % loop through periods 
 	  y = AVAR*y;
 	  DATA(t,:)=y';
end;

time = 1:nimpdat; 
figure(3) 
plot(time, DATA(:,Cpos),'-r');              
hold on; 
plot(time, DATA(:,Kpos), 'b');
plot(time, DATA(:,Zpos),'-k');              
hold off; 
legend('C', 'K', 'Z');

%% Simulation of the model (f)
% First simulate shocks - assume a distribution?
% assume delta~N(0,0.005^2), eps_C~N(0, 0.01^2), eps_Z~N(0,0.01^2)
% realize that all shocks are already in the loglin form!
T=300; %
%T=3000000;
deltat=normrnd(0,sig_D,[1,T]);
epsZt=normrnd(0,sig_Z,[1,T]);
epsCt=normrnd(0,sig_Pr,[1,T]);
shock = zeros(neq,T);
shock([4:6],:)=[deltat;epsZt;epsCt];

y           = zeros(6,1);             % initial value (ie starting at ss)
DATA        = zeros(T,neq);   

for t = 1:T                           % loop through periods 
 	  y = AVAR*y+BVAR*shock(:,t);
 	  DATA(t,:)=y';
end
DATA=DATA(100:T,:);% does not help
C=DATA(:,1);
delta=DATA(:,4);
beta=(delta'*C)/(delta'*delta)
% i.e. higher depreciation (by 100pp) has an effect of -300% on consumption.
% more nicely - 1pp in depr should decrease C by 3%?

% compute Sigma of the variables by hand?
% they are all stationary and mean 0, so Var(yt)=E(yt*yt')=E(yt_1*yt_1')
Epsilon=[sig_D^2,0,0;...
          0, sig_Z^2,0;...
          0, 0, sig_Pr^2];% var cov of the structural shocks
      
Epsilon=[zeros(3,3),zeros(3,3); zeros(3,3), Epsilon];
BEB=BVAR*Epsilon*BVAR';
% W know RHO=A*Sigma, RHO*A'=A*RHO', Sigma=A*RHO'+BEB; Sigma-ASigmaA'=BEB
%Sigma=inv(eye(6)-AVAR)* BEB* inv(AVAR'+eye(6));
% beta_control=Sigma(1,4)/Sigma(4,4)

% Notice that from the structure of AVAR, E(dc) and E(dd) come only from
% the BEB matrix (as there is no intertemporal aspect to delta)
dc=BEB(1,4); 
dd=BEB(4,4);
beta0=dc/dd 

% The estimate is pretty different.
% Which I think makes sense, as C is not iid, in fact the persistance is
% really high, i.e. T has to be really big for C'*d/T to converge to
% E(c*delta).

% If I increase T to be very high, the estimated coefficient is converging
% to beta0.




