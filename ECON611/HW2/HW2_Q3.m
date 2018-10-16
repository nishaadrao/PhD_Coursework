%% ====================================================================== %
% HW2_Q3.m
% The Kalman Filter, Likelihoods, Smoother
%
% A. Yadav 10/15/2018
% ======================================================================= %
addpath ./AIM_subroutines;
clc; clear all; close all;

load('AR_DATA.mat')

%% ======================================================================= %
% (A) Log likelihood for rho=0.8 
% ======================================================================= %


% Parameters
sigma_e = 1;
sigma_m = 0.4;
rho     = 0.8;
v0       = sigma_e/(1-rho^2);

% Run Kalma Filter, store V_x and projections
L = mykalman(0.8,0,v0,XDAT,sigma_m,sigma_e);


%% ======================================================================= %
% (B) Compute MLE estimate of rho 
% ======================================================================= %

rho_vec = 0.75:0.01:0.95;

L_vec   = zeros(length(rho_vec),1);

for p = 1:length(rho_vec)
   
    [L,y_o,y_n, V_o, V_n,V_x] = mykalman(rho_vec(p),0,v0,XDAT,sigma_m,sigma_e);
    
    L_vec(p) = L;
end

% ======================================================================= %
% (C) Kalman smoother
% ======================================================================= %

% Run Kalman recursion for rho=0.75
[L,y_o,y_n, V_o, V_n,V_x] = mykalman(0.99,0,v0,XDAT,sigma_m,sigma_e);

% Create J matrix
J = zeros(length(XDAT)-1,1);

for j=1:length(J)

    J(j) = rho*V_n(j)/V_o(j+1);
    
end

% Create smoothed series
y_smooth      = zeros(length(XDAT),1);
y_smooth(200) = y_n(200);

for j=1:length(J)

    y_smooth(j) = y_n(j) + J(j)*(y_n(j+1)-y_o(j));
    
end

plot(1:200,XDAT,1:200,y_smooth)
legend('x_t','\hat{y}_{t|T}')



function [L,y_o,y_n, V_o, V_n,V_x] = mykalman(rho,y0,v0,X,sigma_m,sigma_e)

    % Initilize some stuff
    y_o      = zeros(length(X),1);
    y_n      = zeros(length(X),1);
    V_o      = zeros(length(X),1);
    V_n      = zeros(length(X),1);
    V_x      = zeros(length(X),1);

    y_o(1)   = y0;
    V_o(1)   = v0;

    % Kalman iteration
    for t = 1:length(X)-1
        
        y_n(t)     = y_o(t) + V_o(t)*(V_o(t) + sigma_m)*(X(t)-y_o(t)); 

        V_x(t)     = V_o(t) + sigma_m;

        V_n(t)     = V_o(t)-V_o(t)^2/V_x(t);
        
        G          = rho*V_o(t)/V_x(t);
        
        y_o(t+1)   = rho*y_o(t) + G*(X(t)-y_o(t));

        V_o(t+1)     = rho^2*V_n(t) + sigma_e;

    end
    
    % Compute probability densities
    f_vec        = zeros(length(y_o),1);

    for t = 1:length(f_vec)

        f_vec(t) = normpdf(y_o(t),V_x(t));

    end

    % Compute log likelihood value
    log_f   = log(f_vec);

    L       = sum(log_f);

end



