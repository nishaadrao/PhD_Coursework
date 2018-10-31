function H = HumanWealth(par, wage, R, tau, dividend, xthis, Hnext)
% Computes the human wealth for grid xthis --  This is just one step of the
% contraction. Call it repeatedly until H converges



global Params;

H = zeros(length(xthis),Params.npp);


%loop over this period
for ip = 1:Params.npp
    
    [~, nthis, sthis] = get_cnbp(xthis,par,R,wage,tau, dividend,ip);
    
    sthis = max(sthis,0);  %small errors around zero cause problems.
    
    H(:,ip) = wage * Params.endow(ip) * nthis;
    
    assets = sthis;
    
    EHprime = 0;
    
    %loop over next period
    for jp=1:Params.npp
        
        pp = transProb(ip,jp);
        if(pp>0)
            
            Hprime = interp1(xthis,Hnext(:,jp),assets,'linear','extrap');
            EHprime = EHprime +  pp * Hprime;
            
            
        end
        
        
    end
    
    H(:,ip) = H(:,ip) + EHprime/R;
    
end
end