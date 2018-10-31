function [update, laborwedge, eulerwedge] = transition_hhld_update(laborwedge, eulerwedge, ii, wage, C, N, ppi, dividend, betashock, parstst, Dstst)

global Params;


eulerwedge0 = eulerwedge;
laborwedge0 = laborwedge;

betapath = Params.beta * betashock;

Rpath = [(1+ii(1:end-1))./ppi(2:end) Params.Rbar];
taupath = Params.B*(1-1./Rpath)/Params.AvgTaxWeight;

parPath = solveback(par2long(parstst),Rpath,wage,taupath,dividend,betapath);
[Cpath, Bpath, Lpath] = simulateforward(Dstst,parPath,Rpath,wage,taupath,dividend);


newlaborwedge = Params.psi1*Lpath.^(Params.psi2)./ (Cpath.^(-Params.sigma).*wage);
neweulerwedge  = Cpath(2:end-2).^(-Params.sigma)./( betapath(2:end-2) .* Rpath(2:end-2) .* Cpath(3:end-1).^(-Params.sigma)) ;




subplot(3,3,1);
plot(ii(1:30)); title('ii')
subplot(3,3,2);
plot(ii(1:10)); title('ii')
subplot(3,3,3);
plot([N(2:end-1)' Lpath(2:end-1)']); title('L')
subplot(3,3,4);
plot([C(2:end-1)' Cpath(2:end-1)']); title('C')
subplot(3,3,5);
plot([eulerwedge(2:end-2)' neweulerwedge']); title('Euler wedge')
subplot(3,3,6);
plot([laborwedge(2:end-1)' newlaborwedge(2:end-1)']); title('labor wedge')
subplot(3,3,7);
plot(Bpath(2:end-1)); title('assets')
drawnow;

RELAX = ones(1,248);
%RELAX(21:end) = 0.98.^([0:227]');  % I used this for ZLB
%RELAX(21:end) = 0.997.^([0:227]');  % I used this for other shocks
laborwedge(2:end-1) = (1-RELAX) .* laborwedge(2:end-1) + RELAX .* newlaborwedge(2:end-1);
eulerwedge(2:end-2) = (1-RELAX(1:end-1)) .* eulerwedge(2:end-2)+  RELAX(1:end-1) .* neweulerwedge;



%Check if we have converged
test = @(orig,new)(max(abs(orig-new)./orig));
update = [test(laborwedge(2:end-1),laborwedge0(2:end-1)) test(eulerwedge(2:end-1), eulerwedge0(2:end-1)) ]
update = max(update);



end



