function X  = DetSimul( X0, eps, fcn, maxit, Dampening, DampeningThresh)
% X  = DetSimul( X0, eps, fcn)
%   Deterministic solution to the model represented by fcn.
%   X0 -- initial guess of solution (nx x T)
%   eps -- exogenous variables (neps x T)
%   fcn -- function that gives residuals of model equations
%   maxit -- maximum number of iterations to try
%   dampening -- initial scale factor of update.
%   DampeningThresh -- residual threshold at which dampening is turned off.
%
%  This program implements the algorithm described in Juillard (1996) 
%  "DYNARE: A program for the resolution and simulation of dynamic models 
%  with forwardd variables through the use of a relaxation algorithm." 
%
% Alisdair McKay
% July 22, 2014


if nargin < 5
    Dampening = 1;
end
if nargin < 6
    DampeningThresh = 0;
end

%initializations
X = X0;

[nx, T] = size(X);
C = zeros(nx,nx,T);
d = zeros(nx,T);
dX = zeros(nx,T);
Fx = zeros(nx,T);





for it = 1:maxit+1
    
    
    
    for t = 2:T-1
        Xt = [X(:,t-1); X(:,t); X(:,t+1)];
        Fx(:,t) = feval(fcn,Xt,eps(:,t));
    end
    
    residual = max(abs(Fx(:)));
    if residual < 1e-8
        disp('done')
        break
    else
        disp(['iteration ' num2str(it) '...residual ' num2str(residual)])
        if residual < DampeningThresh
            Dampening = 1;
        end
            
    end
    
    
    
    
    
    C(:,:,1) = zeros(nx);
    d(:,1) = zeros(nx,1);
    for t = 2:T-1
        Xt = [X(:,t-1); X(:,t); X(:,t+1)];
        [ ~ , S_L, S, S_P] = feval(fcn,Xt,eps(:,t));
%         S_L = adjacob(fcn,Xt,1:nx,nx,eps(:,t));
%         S = adjacob(fcn,Xt,nx+1:2*nx,nx,eps(:,t));
%         S_P = adjacob(fcn,Xt,2*nx+1:3*nx,nx,eps(:,t));
        C(:,:,t) = (S - S_L*C(:,:,t-1))\S_P;
        d(:,t) = -(S - S_L*C(:,:,t-1))\(Fx(:,t) + S_L*d(:,t-1));
    end
    
    dX(:,T) = zeros(nx,1);
    for t = T-1:-1:2
        dX(:,t) = d(:,t) - C(:,:,t)*dX(:,t+1);
    end
    
    
    X = X + Dampening * dX;
    
end

if it > maxit
    warning(['simulation did not converge after ' num2str(maxit) ' iterations.'])
end


end

