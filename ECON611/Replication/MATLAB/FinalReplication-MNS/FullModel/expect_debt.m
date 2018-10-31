function  Edebt  = expect_debt( pvec , pos)
% calculate aggregate debt
%
% pos = true, compute fraction of hhlds with positive assets
% pos = false (default), compute aggregate debt position
%
global Params;
if ~exist('pos','var')
    pos = false;
end
tmp = Params.knotDistrK + Params.borrowCon;
I = tmp>0;
if pos
    tmp = I;
else
    tmp(I) = 0;
    tmp = -tmp;
end
Edebt = expect_k(pvec,1,tmp);


end

