% ECON607_II - HW2 - Q2
%
% Basic RBC model
%
% Anirudh Yadav
% March 18 2018
%
%----------------------------------------------------------------
% 0. Housekeeping (close all graphic windows)
%----------------------------------------------------------------

close all;

%----------------------------------------------------------------
% 1. Defining variables
%----------------------------------------------------------------

var c k i n a y;
varexo e;

parameters alpha beta chi delta rho sigma;

%----------------------------------------------------------------
% 2. Calibration
%----------------------------------------------------------------

load param_nc;
set_param_value('alpha',alpha);
set_param_value('beta',beta);
set_param_value('chi',chi);
set_param_value('delta',delta);
set_param_value('rho',rho);
set_param_value('sigma',sigma);

%----------------------------------------------------------------
% 3. Model
%----------------------------------------------------------------

model;
exp(c)^(-1) = beta*exp(c(+1))^(-1)*(alpha*exp(a(+1))*exp(k)^(alpha-1)*exp(n(+1))^(1-alpha)+1-delta);
0 = exp(c)^(-1)*(1-alpha)*exp(a)*exp(k(-1))^(alpha)*exp(n)^(-alpha)-chi*exp(n);
exp(k) = (1-delta)*exp(k(-1)) + exp(i);
exp(c)+exp(i) = exp(a)*exp(k(-1))^(alpha)*exp(n)^(1-alpha);
a = rho*a(-1) + e;
exp(y) = exp(a)*exp(k(-1))^(alpha)*exp(n)^(1-alpha);
%exp(w) = (1-alpha)*exp(a)*exp(k(-1))^(alpha)*exp(n)^(-alpha);
%exp(r) = alpha*exp(a)*exp(k(-1))^(alpha-1)*exp(n)^(1-alpha);
end;

%----------------------------------------------------------------
% 4. Computation
%----------------------------------------------------------------

initval;
  k = log(10);
  c = log(2.5);
  i = log(0.5);
  a = 0;
  n = log(1/3); 
  e = 0;
end;

shocks;
var e = sigma^2;
end;

steady;

stoch_simul(order = 1,irf=20);





