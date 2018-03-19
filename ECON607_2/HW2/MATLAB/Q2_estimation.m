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

X = [lagmatrix(logtfp,1), t];
[beta,Sigma,e] = mvregress(X,logtfp);



