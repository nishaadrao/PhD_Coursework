function [Cpath, Bpath, Lpath,WealthSharesPath] = simulateforward(D0,parPath,R,wage,tau,dividend)
%[Cpath, Npath] = simulateforward(parPath,ppi,wage,R,labormarketstatus)
%Simulates a distribution of households and computes total consumption,
%labor supply and asset position.
%
% Inputs: D0  -- initial distribution vector
% parPath -- output of solveback
% other -- see solveback for description
%
% Outputs: [X]path -- path for aggregate variable X over transition.


%given today's par and D, compute C, N, B
%find tomorrow's D

global Params;

T = length(R);

Cpath = zeros(1,T);
Bpath = Cpath;
Lpath = Cpath;
WealthSharesPath = zeros(Params.npp,T);
WSP_Summer = sparse(kron(eye(Params.npp),Params.knotDistrK'));

D = D0;

if nargout > 3
    WealthSharesPath(:,1) = WSP_Summer * D;
end


%disp('simulating forward: ' )
for t = 2:T-1
    
    if nargout > 3
        WealthSharesPath(:,t) = WSP_Summer * D;
    end
    
    [ Dprime, Cpath(t), Bpath(t), Lpath(t)] = simulatestep( D,par2wide(parPath(:,t)),R(t),wage(t),tau(t),dividend(t));
    D = Dprime;
    
end
    






