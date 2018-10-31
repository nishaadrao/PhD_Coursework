% File: solveback
% Author: Alisdair McKay
% July 19, 2014
% 
% Description: soles backwards from a period (e.g. steady state) in whcih there
% decision rules are given by parFinal and the prices are wage(end) and so
% on.  R(t) and ppi(t) are the gross nominal interest and inflation from
% date t-1 to t.  wage(t) is the wage in date t.  labormarketstatus(t) is the
% (log) labor market status from t-1 to t.  These four prices should be
% vectors with T elements.  parPath is then n by T where n is the length of
% the output of par2long(par).  wage(end) should be the wage in the steady
% state.
%
% Feel free to use, copy or modify this program in any way.

function parPath = solveback(parFinal,R,wage,tau,dividend,betapath)
%parPath = solveback(parFinal,pricepath)
% 

global Params;

T = length(R);

if ~exist('betapath','var')
    betapath = repmat(Params.beta',1,T);
end


parPath = repmat(par2long(parFinal),1,T);



%disp('solving back: ' )
for t = T-1:-1:2
    
   % if mod(t,50) == 0
   %     disp(num2str(t))
   % end
    
        
    parPath(:,t) = par2long(egm(par2wide(parPath(:,t+1)),R(t),R(t+1),wage(t),wage(t+1),tau(t),tau(t+1),dividend(t),dividend(t+1),betapath(:,t)'));
    
  
    
    
end


end
