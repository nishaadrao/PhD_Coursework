% computes the share of households with wealth <= wstar
function [lwshare, popshare] = expect_low_wealth(pvec,ip,wstar)
    global Params;

    ndk = Params.ndstst;
    I = Params.knotDistrK <= wstar;
        
    offs = (ip-1)*ndk;
    EL = dot( ones(ndk,1), I.* pvec(offs+1:offs+ndk));
    popshare = sum( pvec(offs+1:offs+ndk));
    lwshare = EL/popshare;

    
  
  
      
  
