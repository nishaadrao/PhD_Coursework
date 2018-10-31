function parNew = egm(par,R,Rnext,wage,wagenext,tau,taunext,dividend,dividendnext,beta)
%par = egm_c(par,R,Rnext,wage,wagenext,ppinext,labormarketstatus)
%iterate once on consumption and labor supply using the endog grid method

global Params;


npp = Params.npp;

%% Section 1: solve for savings rule using EGM

% asset grid
xthis = [0; Params.bgrid(2:end)];


nassets = length(xthis);
nc = Params.nc;
assert(nc == nassets);

%compute consumption and marg util for each asset level and income level
MU = NaN(nassets,npp);



for jp=1:npp
    
    
    
    cthis = get_cnbp(xthis,par,Rnext,wagenext,taunext, dividendnext,jp);
    
    
    
    if ~ all(cthis>0)
        disp('negative c')
        jp
        I = cthis <=0;
        find(I)
        cthis(I)
        xthis(I)
        sthis(I)
    end
    assert(all(cthis>0));
    
    MU(:,jp) = margutilC(cthis);
end


%compute expected marg util
MUexp = zeros(nassets,npp);
for ip=1:npp  %loop over previous income states
    for jp = 1:npp %loop over this period income state
        pp = transProb(ip,jp);
        if(pp>0)
            MUexp(:,ip) = MUexp(:,ip) + pp .* MU(:,jp) .* R;
        end
    end
end


clear MU;

Cprev = invmargutilC(bsxfun(@times,beta,MUexp));
assert(all(Cprev(:) > 0))


%at this stage we know the previous c as a function of b'
assert(all(size(Cprev) == [Params.nc, Params.npp]))

%compute labor supply as a function of c and wage
labor = (margutilC(Cprev) .* wage.*repmat(Params.endow,nassets,1)/Params.psi1).^(1/Params.psi2);

%c + b'/R =  b + income + (1-1/R) borrowCon
%b = b'/R + c - income - (1-1/R) borrowCon
bprev = repmat(xthis,1,Params.npp)/R + Cprev - wage.*labor.*repmat(Params.endow,nassets,1)...
    + tau*repmat(Params.tax_weights,nassets,1) - dividend - (1-1/R)*Params.borrowCon;





%pack results back into par
parNew =  NaN(size(par));

for ip = 1:npp

    iConstrained = Params.bgrid <= bprev(1,ip);
    

    
    if any(iConstrained)
        parNew(iConstrained,ip) = egm_solve_constrained(Params.bgrid(iConstrained),ip,wage,tau,dividend,R);
    end
    
    XX = [bprev(:,ip); 1e8];
    tmp = (Cprev(end,ip)-Cprev(end-1,ip))/(bprev(end,ip)-bprev(end-1,ip))*(1e8-bprev(end,ip)) + Cprev(end,ip);
    YY = [Cprev(:,ip); tmp];
    
    parNew(~iConstrained,ip) = interp1(XX,YY,Params.bgrid(~iConstrained),'pchip');
    
   
end



end

function c = egm_solve_constrained(b,ip,wage,tau,dividend,R)

global Params;


n = 0.6 * ones(size(b));
c = b + n*wage*Params.endow(ip) ...
    -tau*Params.tax_weights(ip) + dividend + (1-1/R)*Params.borrowCon;

for it = 1:1000
    f = margutilC(c).*wage.*Params.endow(ip) - Params.psi1 * n.^Params.psi2;
    
    if all(abs(f) < 1e-7)
        break
    end
    
    J = margutilC(c,2)*(wage*Params.endow(ip))^2 - Params.psi1 * Params.psi2 * n.^(Params.psi2-1);
    n = n - f./J;
    c = b + n*wage*Params.endow(ip) + ...
        - tau*Params.tax_weights(ip) + dividend + (1-1/R)*Params.borrowCon;
end

if it == 1000
    error('egm_solve_constrained did not converge')
end

end
