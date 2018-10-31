function  [xe_resid, KL, xaggstst, par, D, Env, X] = compute_stst(xe_in)
%Finds the steady state of the model

    
    global Params modelOptions;

    [margtaxe, taxpaide] = interp_tax(xe_in, Params.incometax);

    
    [ z, M, wage, ii, r, rhat, KL] = SteadyStatePrices( margtaxe );
    
   

    % GIVEN K, SOLVE FOR THE CONSUMPTION FUNCTION BY COLLOCATION:
    % high precision in solution of savings function_ is required to achieve
    %   converge in outer loop!
    [par,check] = broydn(@eulerres_stst,Params.parstart,[1e-11,0,1],1+ii,wage);
    if(check~=0)
        %save broyError par;
        warning('compute_stst:broyerror','broydn not converged');
    end
    
    
    Params.parstart = par;
    par = par2wide(par);

    labmktstat = 0;
    
    Pi = forwardmat(1,par,labmktstat);
    D = invdistr(Pi);
    % D = invd2(Pi);
    clear Pi;





    xaggstst = agg_stst(margtaxe, taxpaide, D,par);

    xe_resid = xe_in - xaggstst(strcmp(Params.aggnames,'xe'));
    
    if(nargout>4)
        Env.Dstst = D;
        Env.R = 1+r;
        Env.wage = wage;
        Env.parstst = par;
        Env.aggstst = xaggstst;
    end
    
    if nargout > 6
        
        if modelOptions.search
            X = [zeros(Params.nStatesDistr-1,1);xaggstst;par2long(par)];
        else
            parshort = par2long(par);
            parshort = parshort(Params.npp*Params.nv+1:end);  %drop the value function parameters
            X = [zeros(Params.nStatesDistr-1,1);xaggstst;parshort];
        end
  
    end



end

