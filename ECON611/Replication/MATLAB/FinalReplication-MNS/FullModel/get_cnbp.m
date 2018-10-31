function [c, n, bp] = get_cnbp(x,par,R,wage,tau, dividend,incomeIndex)

global Params;

c = interp1(Params.bgrid,par(:,incomeIndex),x,'pchip');

if nargout > 1
n = (margutilC(c)*wage*Params.endow(incomeIndex)/Params.psi1).^(1/Params.psi2);
bp = R*(x + wage*n*Params.endow(incomeIndex) - tau*Params.tax_weights(incomeIndex) + dividend + (1-1/R)*Params.borrowCon - c);
end
