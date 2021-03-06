%% ECON607_II - HW2 - Q2
%  
%  Initialization file for Dynare file 'Q2.mod', which implements basic RBC
%  model
% 
%  Anirudh Yadav
%  18 March 2018
%
%  NOTES:...
%
%----------------------------------------------------------------
% 0. Housekeeping (close all graphic windows)
%----------------------------------------------------------------

clear;
close all;

%----------------------------------------------------------------
% 1. Parameter values
%----------------------------------------------------------------

alpha = 1/3;
beta = 0.99;
delta = 0.025;
v = 0.72;
%v = 0.72;
rho = 0.97; % check with others
sigma = 0.0085; % check with others!

%----------------------------------------------------------------
% 2. Compute chi to target N*=1/3
%----------------------------------------------------------------

% Compute stead-state capital-labour ratio
knratio = (alpha/(1/beta - (1-delta)))^(1/(1-alpha));

% Compute chi to target N*=1/3
nbar = 1/3;
chi = (((1-alpha)*knratio^alpha)/(knratio^alpha - delta*knratio))/(nbar^(1+1/v));

%----------------------------------------------------------------
% 3. Save parameter values for Dynare
%----------------------------------------------------------------

save param_nc alpha beta chi rho sigma v

%----------------------------------------------------------------
% 4. Run Dynare
%----------------------------------------------------------------

dynare Q2

%----------------------------------------------------------------
% 5. Manually compute moments to check if they match Dynare output tables
%----------------------------------------------------------------

% With v=0.72, dynare spits out std(y) = 0.013514 and std(c) = 0.005356
% Let's see if we can match these with manual calculations using the
% simulated data.

% y check
ytilde = y - hpfilter(y,1600);
std(ytilde)

% c check
ctilde = c - hpfilter(c,1600);
std(ctilde)

% They're pretty close (std ytilde is a bit higher than Dynare's output).
