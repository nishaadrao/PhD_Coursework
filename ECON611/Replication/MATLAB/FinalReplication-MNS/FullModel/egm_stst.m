function par = egm_stst(par,R,wage,tau,dividend)
  %compute steady state savings and labor supply using endog grid method
  
%  global Params;
  
T = 30;  
maxit = 500;
for it = 1:maxit
    parPath = solveback(par2long(par),repmat(R,T,1),repmat(wage,T,1),repmat(tau,T,1),repmat(dividend,T,1));
    parPath = parPath(:,2:end);
    [test, itest] = max(abs(parPath(:,1)-parPath(:,2)));
    disp(['Update to policy rule: ' num2str(test)]);
    %disp(['itest: ' num2str(itest)])
    par = par2wide(parPath(:,1));
    if test < 1e-7, break; end
    
%     parold = par2wide(parPath(:,2));
%     parold3 = par2wide(parPath(:,3));
%     plot(Params.bgrid,[par(:,1) parold(:,1) parold3(:,1)]);
%     drawnow;
%     %pause
end



    
    
end

