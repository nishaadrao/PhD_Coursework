% ================================================================================= % 
% irbc_N_bond.m 
%
% Modification of irbc_new.m to restrict assets to only uncontingent bonds. 
%
% C. House 9/26/2018
% ================================================================================= % 

addpath ./AIM_subroutines;
clear all;
clc;
close all; 

N = 2;              % Number of countries.

nimp_a 	= 100;		% number of years in impulse response
nplot_a = 100;    	% number of years to be plotted
h       = 1/12;     % fraction of a year per period (h=1 is one year)
year_plot = 5;      % how many years to plot

N_plot  = [1 2];    % which countries to plot

nimp    = nimp_a/h;   
nplot   = nplot_a/h;  

% Model Parameters 
a       = .35; 		% capital share
s       = 1;		% intertemporal elasticity of substitution (s=1:log u)
eta     = 1;        % Frisch labor supply elasticity
zeta    = 1/h;      % investment adjustment cost
% zeta    = 0; 
d_ann   = .10; 		% annual depreciation rate
r_ann   = .04;      % annaul discount rate

beta    = exp(-r_ann*h);        % per period discount factor. 
d       = 1 - exp(-h*d_ann);    % per period depreciation rate

% Technology Process

%rho_ann 	=  .05^12;  % annual persistence
rho_ann 	=  .5;      % annual persistence

rho 	=  rho_ann^h;   % per period persistence
rho1 	=  0;           % per period spillover from i to i+1
rho2 	=  0;           % per period spillover from i to i-1

GAMMA = zeros(N,N);

GAMMA(1,1) = rho;
GAMMA(1,2) = rho1;
GAMMA(1,N) = rho2;

GAMMA(N,N) = rho;
GAMMA(N,1) = rho1;
GAMMA(N,N-1) = rho2;

for i = 2:N-1;
  GAMMA(i,i-1:i+1) = [rho2 rho rho1];
end;

if max(abs(eig(GAMMA))) >= 1;
  disp('unstable technology process');
  return;
end;

asigma = .01;

shockvar = eye(N)*asigma^2;	% uncorrelated tech shocks.

% ======================================================================= %
% Nonstochastic steady state:                                    
% ======================================================================= %

K   = (a/((1/beta) - 1 + d))^(1/(1-a));
Y   = K^a; 
I   = d*K; 
C   = Y - I; 
W   = (1-a)*Y;
R   = a*Y/K;

%******************************************************************
%**** Begin setup for constructing log-linear system coef matrix **
%******************************************************************

nlead = 1;  % Number of leads in system 
nlag  = 1;  % Number of lags in system 

xnum = (11*N)+1; 	%  Number of variables in system 
neq = xnum;		%  Number of equations (same)

Cpos 	= [];
Lpos 	= [];
Kpos 	= [];
Ipos    = []; 
Ypos 	= [];
Zpos 	= [];
e_Zpos 	= [];

for i = 1:N;
  
  Cpos(i)   =  (10*(i-1)) + 1;
  Kpos(i)   =  (10*(i-1)) + 2;
  Lpos(i)   =  (10*(i-1)) + 3;
  Ypos(i)   =  (10*(i-1)) + 4;
  Ipos(i)   =  (10*(i-1)) + 5;
  qpos(i)   =  (10*(i-1)) + 6;
  Wpos(i)   =  (10*(i-1)) + 7;
  Rpos(i)   =  (10*(i-1)) + 8;
  Zpos(i)   =  (10*(i-1)) + 9;
  Spos(i)   =  (10*(i-1)) + 10; 
 
end; 

% shadow value of consumption.

rpos = (N*10) + 1; 

% there are N structural shocks -- ordered last.

for i = 1:N;
  e_Zpos(i)  = rpos + i; 
end; 

colzero = 0+nlag*xnum;      % Position counter for start of contemp. coefs 
                            % i.e. if there is one lag and 5 equations then the "first"
                            % contemporaneous variable is at position 6 (6 = 1*5+cpos)
collead = 0+nlag*xnum+xnum; % Position counter for start of lead coefs 
collag  = 0;                % Position counter for start of lag coefs  

% Indicators for contemporanous coefs for each variable: 

Czero       = colzero+Cpos;
Kzero       = colzero+Kpos;
Lzero       = colzero+Lpos;
Yzero       = colzero+Ypos;
Izero       = colzero+Ipos;
qzero       = colzero+qpos;
Wzero       = colzero+Wpos;
Rzero       = colzero+Rpos;
Zzero       = colzero+Zpos;
Szero       = colzero+Spos;
rzero       = colzero+rpos; 
e_Zzero     = colzero+e_Zpos;

% Indicators for lead coefficients for each variable: 

Clead       = collead+Cpos;
Klead       = collead+Kpos;
Llead       = collead+Lpos;
Ylead       = collead+Ypos;
Ilead       = collead+Ipos;
qlead       = collead+qpos;
Wlead       = collead+Wpos;
Rlead       = collead+Rpos;
Zlead       = collead+Zpos;
Slead       = collead+Spos;
rlead       = collead+rpos; 
e_Zlead     = collead+e_Zpos;

% Indicators for lag coefficients for each variable: 

Clag       = collag+Cpos;
Klag       = collag+Kpos;
Llag       = collag+Lpos;
Ylag       = collag+Ypos;
Ilag       = collag+Ipos;
qlag       = collag+qpos;
Wlag       = collag+Wpos;
Rlag       = collag+Rpos;
Zlag       = collag+Zpos;
Slag       = collag+Spos;
rlag       = collag+rpos; 
e_Zlag     = collag+e_Zpos;

% NOTE: these are all vectors of position counters
% Now I have one vector with all of the leads, lags, etc in one column

% Determine number of coefficients per equation: 

ncoef = neq*(nlag+nlead+1);

cof = zeros(neq,ncoef);             % Coef matrix --- Each row is an equation

% ================================================================================= % 
% ================================================================================= % 
% Fill in the coefficient matrix (cof)
% ================================================================================= % 
% ================================================================================= % 

% ================================================================================= % 
% Euler equations
% ================================================================================= % 

for i = 1:N;

  cof(i,Czero(i))   = 1/s;
  cof(i,Clead(i))   = -1/s;
  cof(i,rzero)      = beta;
  
end;

% ================================================================================= % 
% Labor supply
% ================================================================================= % 

for i = 1:N;

  cof(N+i,Lzero(i)) = -(1/eta);
  cof(N+i,Wzero(i)) = 1;
  cof(N+i,Czero(i))  = -(1/s);
  
end;

% ================================================================================= % 
% Investment foc's
% ================================================================================= % 

for i = 1:N;

  cof((2*N)+i,qzero(i)) = -1;
  cof((2*N)+i,Czero(i))  = -1/s;
  cof((2*N)+i,Izero(i)) = zeta*d;
  cof((2*N)+i,Klag(i))  = -zeta*d;

end;

% ================================================================================= % 
% Wages
% ================================================================================= % 

for i = 1:N; 
    
    cof((3*N)+i,Wzero(i))   = -1;
    cof((3*N)+i,Zzero(i))   = 1;
    cof((3*N)+i,Klag(i))    = a;
    cof((3*N)+i,Lzero(i))   = -a;

end;

% ================================================================================= % 
% Rental Rate on Capital
% ================================================================================= % 

for i = 1:N; 
    
    cof((4*N)+i,Rzero(i))   = -1;
    cof((4*N)+i,Zzero(i))   = 1;
    cof((4*N)+i,Klag(i))    = a-1;
    cof((4*N)+i,Lzero(i))   = 1-a;

end;

% ================================================================================= % 
% Shadow value of Capital
% ================================================================================= % 

for i = 1:N; 
    
    cof((5*N)+i,qzero(i))   = -1;
    cof((5*N)+i,Rlead(i))   = beta*R;
    cof((5*N)+i,Clead(i))   = -beta*R/s;
    cof((5*N)+i,Kzero(i))   = -beta*zeta*(d^2);
    cof((5*N)+i,Ilead(i))   = beta*zeta*(d^2);
    cof((5*N)+i,qlead(i))   = beta*(1-d);

end;


% ================================================================================= % 
% Production
% ================================================================================= % 

for i = 1:N; 
    
    cof((6*N)+i,Yzero(i))   = -1;
    cof((6*N)+i,Zzero(i))   = 1;
    cof((6*N)+i,Klag(i))    = a;
    cof((6*N)+i,Lzero(i))   = 1-a;

end;


% ================================================================================= % 
% Capital Accumulation
% ================================================================================= % 

for i = 1:N; 
    
    cof((7*N)+i,Kzero(i))   = -1;
    cof((7*N)+i,Izero(i))   = d;
    cof((7*N)+i,Klag(i))    = 1-d;

end;


% ================================================================================= % 
% Resource Constraint
% ================================================================================= % 

for i = 1:N; 
    
    cof((8*N)+i,Yzero(i))   = -1;
    cof((8*N)+i,Slag(i))    = -1/beta;
    cof((8*N)+i,Izero(i))   = I/Y;
    cof((8*N)+i,Czero(i))   = C/Y;
    cof((8*N)+i,Szero(i))   = 1;

end;

count = 9*N;

% ================================================================================= % 
% Bond Market Clearing
% ================================================================================= % 

for i = 1:N; 
    
    cof(count+1,Szero(i))   = 1;

end;

count = count + 1;


% ================================================================================= % 
% Shock Processes
% ================================================================================= % 

cof(count+1,Zzero(1))   = 1;
cof(count+1,Zlag(1))    = -rho;
cof(count+1,Zlag(N))    = -rho1;
cof(count+1,Zlag(2))    = -rho2;
cof(count+1,e_Zzero(1)) = -1;

for i = 2:N-1;
    
    cof(count+i,Zzero(i))   = 1;
    cof(count+i,Zlag(i))    = -rho;
    cof(count+i,Zlag(i-1))    = -rho1;
    cof(count+i,Zlag(i+1))    = -rho2;
    cof(count+i,e_Zzero(i)) = -1;

end;

cof(count+N,Zzero(N))   = 1;
cof(count+N,Zlag(N))    = -rho;
cof(count+N,Zlag(N-1))    = -rho1;
cof(count+N,Zlag(1))    = -rho2;
cof(count+N,e_Zzero(N)) = -1;


count = count+N; 

% ================================================================================= % 
% Structural Shocks. 
% ================================================================================= % 

for i = 1:N;
  cof(count+i,e_Zzero(i)) = 1;
end;


%*****************************************************************
%*************Begin Solution Algorithm ***************************
%*****************************************************************
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
%    nlag     Number of lags.                                      
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
%*****************************************************************

% Use AIM procedure to solve model: 
uprbnd = 1+1e-8;    % Tolerance values for AIM program 
condn = 1e-8;

% ================================================================================= % 
% Run AIM
% ================================================================================= % 

[cofb,rts,ia,nex,nnum,lgrts,mcode] = ...
       aim_eig(cof,neq,nlag,nlead,condn,uprbnd);

scof = obstruct(cof,cofb,neq,nlag,nlead);

% need to calculate amat and b
% ===============================
s0 = scof(:,(neq*nlag+1):neq*(nlag+1)); 	%Contemp. coefs from obs. structure
amat=zeros(neq*nlag,neq*nlag);   		% Initialize A matrix 
bmat=cofb(1:neq,((nlag-1)*neq+1):nlag*neq);  	% Lag 1 coefficients 
i=2;
while i<=nlag;
  bmat=[bmat cofb(1:neq,((nlag-i)*neq+1):(nlag-i+1)*neq)]; % Lag i coefs 
  i=i+1;
end;
amat(1:neq,:)=bmat;  				% Coefs for equations 
if nlag>1;
 amat((length(cofb(:,1))+1):length(amat(:,1)),1:neq*(nlag-1))=eye(neq*(nlag-1));
end;
b = zeros(length(amat(:,1)),length(s0(1,:)));
b(1:length(s0(:,1)),1:length(s0(1,:))) = inv(s0);  % Store coefs 


% ================================================================================= % 
% VAR form matricies. 

AVAR = amat; 
BVAR = b; 

% ================================================================================= % 

% ================================================================================= % 
% Impulse Response Function
%
% Technology shock to country #1
% ================================================================================= % 


shock = zeros(neq,1);                           % Shock vector 
shock(e_Zpos(1),1) = 1;                         % Shock variable, size of shock 
impdat = impf(AVAR,BVAR,shock,nimp,neq);          % Call impulse respons proc 
dat=(1:nplot)';                                 % Date variable for plotting 

M = length(N_plot); 

yr      = 1/h;      % number of "clicks" in a year. 
yr_count = 0:yr*year_plot:nplot; 
yr_label = yr_count*h; 

colordef white
orient landscape

for j = 1:length(N_plot); 

 figure(j)
 plot(dat,impdat(1:nplot,Ypos(N_plot(j))),'g-', 'LineWidth', 1.5); 
 hold on; 
 plot(dat,impdat(1:nplot,Ipos(N_plot(j))),'k-', 'LineWidth', 1.5);
 plot(dat,impdat(1:nplot,Lpos(N_plot(j))),'r-', 'LineWidth', 1.5);
 plot(dat,impdat(1:nplot,Cpos(N_plot(j))),'m-', 'LineWidth', 1.5); 
 hold off; 
 legend('GDP','I','L','C','Location','NorthEast');
 set(gca,'Xtick',yr_count);
 set(gca,'Xticklabel',yr_label);
 set(gca,'YGrid','on')

end; 


 colordef white
 figure(M+1) 
 plot(dat,impdat(1:nplot,Cpos(N_plot(1))),'g-'); 
 hold on; 
 for j = 2:M; 
     plot(dat,impdat(1:nplot,Cpos(N_plot(j))),'b-');
 end; 
 hold off; 
 title('Consumption')
 set(gca,'Xtick',yr_count);
 set(gca,'Xticklabel',yr_label);
 set(gca,'YGrid','on')


 colordef white
 figure(M+2) 
 plot(dat,impdat(1:nplot,Spos(N_plot(1))),'g-'); 
 hold on; 
 for j = 2:M; 
     plot(dat,impdat(1:nplot,Spos(N_plot(j))),'b-');
 end; 
 hold off; 
 title('Savings')
 set(gca,'Xtick',yr_count);
 set(gca,'Xticklabel',yr_label);
 set(gca,'YGrid','on')


 