%% ECON607_II - HW3 - Q1(4)
%
%  Anirudh Yadav
%  29 March 2018
%
%
%----------------------------------------------------------------
% 0. Housekeeping (close all graphic windows)
%----------------------------------------------------------------

clear;
close all;


%----------------------------------------------------------------
% 1. Calibration
%----------------------------------------------------------------
beta = 0.996;
tau = 0.4;
chi = 0.034;
phi = 0.5;
gamma = 0.399;

syms theta c n mu;

f = 2.32*theta^(1/2);

mu = f/theta;

c = (f-chi)/(chi+f);

w = phi*(1+theta) + (1-phi)/(1-tau)*gamma*c;

eqn1 = (mu + 1 - chi -mu*w)*beta == 1;

eqn2 = (theta*(1-beta*(1-chi)))/(beta*f) - 1 + phi*(1+theta) + ((1-phi)*gamma*(f-chi*theta))/((1-tau)*(f+chi))== 0;

soltheta1 = solve(eqn1,theta);
soltheta2 = solve(eqn2,theta);


% Get solution for theta
theta1 = vpa(soltheta1);
theta2 = vpa(soltheta2);

% Drop symbols
clear c n mu;

% Compute other steady state values
f = 2.32*theta1^(1/2);

mu = f/theta1;

c = (f-chi)/(chi+f);

w = phi*(1+theta1) + ((1-phi)/(1-tau))*gamma*c;

n = f/(chi+f);

v = chi/mu;


