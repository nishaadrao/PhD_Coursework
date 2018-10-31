function  parstst = steadystatepolicy(Y)


%STARTING VALUES FOR_ CONSUMPTION FUNCTION AND LABOR SUPPLY:


global Params;


if isfield(Params,'parstart')
    par0 = Params.parstart;
else
    %initial guess
    pars = 0.3+[0.1*Params.bgrid]; %first element is sup {a :a'(a) = 0}


    par0 = repmat(pars,1,Params.npp);

end



if size(par0,2) == 1
    par0 = par2wide(par0);
end

%to get a good initial guess, use endog grid method but it does not need to
%converge fully
%disp('Using endog. grid method to find policy rules.');

[R, w, tau,dividend] = steadystateprices(Y);
parstst = egm_stst(par0,R,w,tau,dividend);

Params.parstart = parstst;




end

