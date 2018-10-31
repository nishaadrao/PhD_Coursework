
% Computes the effective labor supply,
% Inputs:
%   pvec:  vector of probabilities (histogram weights)
%   par: policy rule parameters
%   prodWeight: if true (default), weight hours worked by productivity.

function EL = expect_L(pvec,par,prodWeight,R, wage,tau,dividend)
  global Params;
  
  ndk = Params.ndstst;
  
  EL = 0;
  
  for ip=1:Params.npp
      
      
      if  prodWeight
          prod_scale = Params.endow(ip);
      else
          prod_scale = Params.endow(ip)>1e-6;
      end
     
      
      [~, nsup] = get_cnbp(Params.knotDistrK,par,R,wage,tau, dividend,ip);
      
      offs = (ip-1)*ndk;
      EL = EL + sum(pvec(offs+1:offs+ndk) .* nsup * prod_scale);
  end
  
  assert(ip*ndk == length(pvec));
      
  
