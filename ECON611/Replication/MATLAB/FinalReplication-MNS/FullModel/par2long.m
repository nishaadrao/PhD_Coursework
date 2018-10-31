function par = par2long(par)



global Params;


nc = Params.nc;
npp = Params.npp;


[m n] = size(par);


if (m == nc * npp && n == 1)
    return
end

assert(m == nc && n == npp);

par = par(:);



