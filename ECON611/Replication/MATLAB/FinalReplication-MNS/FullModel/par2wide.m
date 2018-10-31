function par = par2wide(par)

    

    global Params;


    nc = Params.nc;
    
    npp = Params.npp;

    
    [m n] = size(par);
    
    
        if m == nc && n == npp
            return
        end
    

if (m == nc* npp && n == 1)

    par = reshape(par,nc,npp);

else
    error(['format of par supplied to par2wide is not recognized' num2str(m) ' ' num2str(n)])
end

