% [EC Emargtax Etaxpaid] = expect_C(pvec,par,R,wage)
% Computes the aggregate consumption and taxes of households
% Inputs:
%   pvec:  vector of probabilities (histogram weights)
%   par: policy rule parameters
%   R: interest rate
%   wage
%   disc: vector of discrete states levels over which to integrate
%
% Alisdair McKay 2/24/12

function [EC, EL] = expect_C(pvec,par,R,wage,tau,dividend)

  global Params;

  
  
  npp = Params.npp;
  par = par2wide(par);
  
  ndk = Params.ndstst;
  
  EC = 0;
  EL = 0;

  
  dgrid = 1:npp;
  
    
  %loop over discrete grid
  for ip=dgrid
      
      [c, n] = get_cnbp(Params.knotDistrK,par,R,wage,tau, dividend,ip);
      
      
      
      
      %build up expectation
      offs = (ip-1)*ndk;
           
      EC = EC + sum(pvec(offs+1:offs+ndk) .* c);
      EL = EL + sum(pvec(offs+1:offs+ndk) .* n .* Params.endow(ip));
      
      
  end
      
  
end
