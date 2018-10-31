function eqm  = SolveForTransition( Rpath ,  wagepath, dividendpath, Spath, steadystate, parstst, Dstst)
%Iteratively solves for transition path for section 3

global Params;

T = length(Rpath);

taupath = Params.B*(1-1./Rpath)/Params.AvgTaxWeight;

fwage = figure;

for outer_it = 1:100
    
    for inner_it = 1:100
        
        
        
        parPath = solveback(par2long(parstst),Rpath,wagepath,taupath,dividendpath);
        [Cpath, ~, Lpath,WealthSharesPath] = simulateforward(Dstst,parPath,Rpath,wagepath,taupath,dividendpath);
        Y = Cpath;
        N = Spath.*Y;
        subplot(2,1,1);
        plot([Lpath(2:T-1)'/steadystate.Y N(2:T-1)'/steadystate.Y ones(size(N(2:T-1)'))])
        subplot(2,1,2);
        plot([wagepath' steadystate.w*ones(size(wagepath'))])
        drawnow
        
        
        %adjust wage, dividend
        oldwage = wagepath;
        wagepath(2:T-1) = wagepath(2:T-1).*(N(2:T-1)./Lpath(2:T-1)).^Params.psi2;
        
        %wagepath = 0.95 * wagepath + 0.05 * oldwage;
        wagepath = 0.25 * wagepath + 0.75 * oldwage;
        
        dividendpath(2:end-1) = (Y(2:T-1) - wagepath(2:T-1).*N(2:T-1));
        
        test = max(abs(wagepath./oldwage - 1));
        disp(['Residual in wage path: ' num2str(test)])
        if test < 1e-6, break; end
        
        
    end
    
    %Now compute inflation and S
    %solve backwards from the end for pi
    %then solve forwards from the beginning for S
    
    % initialize these with steady state values
    pbarA = Params.mu * steadystate.w * steadystate.Y /  (1-Params.beta(end)*(1-Params.theta));
    pbarB = steadystate.Y/  (1-Params.beta(end)*(1-Params.theta));
    
    %solve for inflation
    ppipath = ones(1,T);
    for t = T-1:-1:2
        pbarA = Params.mu*wagepath(1+t) * Y(t) + Params.beta(end)*(1-Params.theta)* ppipath(t+1)^(- Params.mu/(1- Params.mu))*pbarA;
        pbarB = Y(t) + Params.beta(end)*(1-Params.theta) * ppipath(t+1)^(- 1/(1- Params.mu))*pbarB;
        pstar(t) = pbarA/pbarB;
        ppipath(t) = ((1-Params.theta)/(1-Params.theta*pstar(t)^(1/(1-Params.mu))))^(1-Params.mu);
    end
    
    % solve for S
    oldS = Spath;
    Spath = ones(1,T);
    Slast = 1;
    for t = 2:T-1
        Spath(t) = (1-Params.theta)*Slast*ppipath(t)^(-Params.mu/(1-Params.mu)) + Params.theta*pstar(t)^(Params.mu/(1-Params.mu));
        Slast = Spath(t);
    end
    
    test = max(abs(Spath./oldS - 1 ));
    disp(['Residual in S path: ' num2str(test)])
    if test < 1e-6, break; end
end

eqm.S = Spath;
eqm.w = wagepath;
eqm.ppi =ppipath;
eqm.Y = Y;
eqm.R = Rpath;
eqm.tau = taupath;
eqm.dividend = dividendpath;
eqm.WealthShares = WealthSharesPath;

end

