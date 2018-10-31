function [ resid ] = check_steady_state( x )

global Params;

if Params.beta_heterog
    Params.beta = kron(x(1)/100+[-Params.betad 0],ones(1,Params.npp/Params.nbeta));
else
    Params.beta = x(1)/100;
end
Y = x(2);


[R, w, tau, dividend] = steadystateprices(Y);


parstst = steadystatepolicy(Y);
% find stationary distribution of wealth
Pi = forwardmat(1,parstst,R,w,tau,dividend);
D = invdistr(Pi);
clear Pi;
aggassets = expect_k(D);

[EC, EL] = expect_C(D,parstst,R,w,tau,dividend);



%the first target is assets over GDP
% the second one is the aggregate resource constraint
resid = [aggassets/EC - Params.asset_target;  EC-EL]
%error('stop here')

% 
% 
% b = Params.knotDistrK + Params.borrowCon;
% b = repmat(b,Params.npp,1);
% [~, Lorenz] = gini(D, b, false);
% Lorenz80 =  1 - Lorenz(find(Lorenz(:,1) <= 0.8,1,'last'),2);
% 
% 
% 
% %the first target is assets over GDP
% resid = [aggassets/EC - Params.asset_target;  Lorenz80 - Params.LorenzTarget; EC-EL]


end

