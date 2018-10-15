% ======================================================================= %
% HW_2_Q2.m
% A 2-agent endowment economy with non-contingent bond
%
% A. Yadav 10/15/2018
% ======================================================================= %
addpath ./AIM_subroutines;
clc; clear all; close all;

N = 2;              % Number of agents

nimpdat 	= 200;	% No. of quarters in impulse response. 

% ======================================================================= %
% 			     Parameters 
% ======================================================================= %

% annual rates
r           = .02;		% quarterly discount rate (4 pct annual)
beta        = 1/(1+r); 


% Shock Variances
sig_Y       = .01; 


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

xnum    = 8;
neq 	= xnum;

% ======================================================================= %
% Position Counters 
% ======================================================================= %

Cpos 	= [];
Ypos 	= [];
Spos 	= [];
e_Ypos 	= [];

for i = 1:N;
  
  Cpos(i)   =  (3*(i-1)) + 1;
  Ypos(i)   =  (3*(i-1)) + 2;
  Spos(i)   =  (3*(i-1)) + 3;
  e_Ypos(i) =  (3*(i-1)) + 4;
 
end; 


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
Szero       = colcount + Spos;
e_Yzero	    = colcount + e_Ypos;

% =========================================================
% Indicators for lead coefficients for each variable: 

colcount  = collead;

Clead		= colcount + Cpos; 
Ylead		= colcount + Ypos;
Slead       = colcount + Spos;
e_Ylead	    = colcount + e_Ypos;

% =========================================================
% Indicators for lag coefficients for each variable: 

colcount  = collag;

Clag		= colcount + Cpos; 
Ylag		= colcount + Ypos;
Slag        = colcount + Spos;
e_Ylag	    = colcount + e_Ypos;

% ==================================================================================== %
% Coefficient Matrix
% ==================================================================================== %

% Number of coefficients per equation: 
ncoef   = neq*(nlag+nlead+1);

cof     = zeros(neq,ncoef);         % Coef matrix. Each row is an equation 

% Euler Equations

for i = 1:N;

  cof(i,Czero(i))   = 1;
  cof(i,Clead(i))   = -1;
  
end;

% Budget constraints

for i = 1:N; 
    
    cof(N+i,Czero(i))   = -1;
    cof(N+i,Szero(i))   = -1;
    cof(N+i,Yzero(i))   = 1;
    cof(N+i,Slag(i))    = 1+r;

end;

% Income processes

for i = 1:N; 
    
    cof(2*N+i,Yzero(i))   = -1;
    cof(2*N+i,e_Yzero(i))   = 1;

end;

% Bond market clearing

for i = 1:N; 
    
    cof(6,Szero(i))   = 1;

end;


% SHOCK: Y

for i = 1:N;
  cof(6+i,e_Yzero(i)) = 1;
end;

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

%}
%% ====================================================================== %
% ======================================================================= %
%                         IMPULSE RESPONSES
% ======================================================================= %
% ======================================================================= %
%{
% Income shock to agent 1

shock              = zeros(neq,1);              % Shock vector 
shock(e_Ypos(1),1) = 1;                         % Shock variable, size of shock 
time               =(1:nimpdat)';               % Date variable for plotting


impdat = impf(AVAR,BVAR,shock,nimpdat,neq);        % Call impulse respons proc 

figure(1)
plot(time,impdat(1:nimpdat,Cpos(1)),'g-'); 
hold on;
plot(time,impdat(1:nimpdat,Cpos(2)),'b-');
hold off;
title('Consumption')
legend('c^A','c^B');

figure(2)
plot(time,impdat(1:nimpdat,rpos));
%}




