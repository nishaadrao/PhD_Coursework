function [ Dprime, C, Assets, L ] = simulatestep( D,par,R,wage,tau,dividend)
%[ Dprime, C, Assets ] = simulatestep( D,par,ppi,wage,R,labormarketstatus)
%
% Inputs;
% D is distribution of assets entering date t before adjusting for inflation
% par is labor supply in t and savings from t to t+1
% ppi is inflation from t-1 to t
% wage is wage in t
% R is gross nominal interest from t-1 to t
% labormarketstatus describes transitions from t to t+1
%
% Outputs:
% Dprime is distribution entering t+1 before adjusting for inflation
% C is agg impatient cons in t
% Assets is agg impatient saving from t to t+1
% tax payment is agg impatient income tax payment in t




%household consumption
[C, L] =  expect_C(D,par,R,wage,tau,dividend);

Pi = forwardmat(1,par,R,wage,tau,dividend);
Dprime = forward(D,Pi);


Assets = expect_k(Dprime);  % def of assetsh

assert(abs(sum(Dprime)-1) < 1e-6)

end

