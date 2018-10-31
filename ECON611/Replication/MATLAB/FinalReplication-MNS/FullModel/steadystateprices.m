function [R, w, tau, dividend] = steadystateprices(Y)

global Params;

R = Params.Rbar;
w = 1/Params.mu;
if nargout > 2
    tau = Params.asset_target * Y *(1-1/R)/Params.AvgTaxWeight;
    dividend = Y*(1-w);
end