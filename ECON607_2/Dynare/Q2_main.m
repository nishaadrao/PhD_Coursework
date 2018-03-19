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
rho = 0.95; % to update!
sigma = 1; % to update!

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

save param_nc alpha beta chi rho sigma

%----------------------------------------------------------------
% 4. Run Dynare
%----------------------------------------------------------------

dynare Q2
