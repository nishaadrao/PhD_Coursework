%% ECON607_II - HW2 - Q2 - first part
%  
%  File for estimating persistence of productivity process (rho) and
%  variance of productivty shocks (sigma)
% 
%  Anirudh Yadav
%  19 March 2018
%
%  NOTES:...
%
%----------------------------------------------------------------
% 0. Housekeeping (close all graphic windows)
%----------------------------------------------------------------

clear;
close all;

%----------------------------------------------------------------
% 1. Import and prepare data
%----------------------------------------------------------------

% Import raw data
  
  tfp = xlsread('/Users/Anirudh/Desktop/PhD/ECON607/HW2/Data/tfp.xlsx');

% Log the data

  logtfp = log(tfp(:,2));

% Create a time variable

  t = (1:1:length(logtfp))';
  
%----------------------------------------------------------------
% 2. Run OLS with linear time trend
%----------------------------------------------------------------

% Create a vector of ones
intercept = ones(280,1);

% Create data matrix
X_1 = [intercept, t];

% Run OLS, store estimates
[rho_a,Sigma_A,e] = mvregress(X_1,logtfp);

%----------------------------------------------------------------
% 3. Run AR(1) on residuals
%----------------------------------------------------------------

% Create data matrix
X_2 = [lagmatrix(e,1), intercept];

% Run OLS, store estimates
[rho_e,Sigma_e,u] = mvregress(X_2,e);

%----------------------------------------------------------------
% 4. Compute std dev of residuals from AR(1) regression
%----------------------------------------------------------------
 
% Exclude first entry which is NaN
sigma=std(u(2:end));