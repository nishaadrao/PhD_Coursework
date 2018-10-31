close all
clear all
clc

setpath();

run('FullModel/parameters/initparams');


%% Adjust the discount factor to target a certain level of assets


Y0 = 0.6;
optset('broyden','tol',1e-5);
x = broyden(@check_steady_state,[100*Params.beta(end);Y0]);
if Params.beta_heterog
    Params.beta = kron(x(1)/100+[-Params.betad 0],ones(1,Params.npp/2));
else
    Params.beta = x(1)/100;
end
steadystate.Y = x(2);
Params.B = steadystate.Y * Params.asset_target;

assert(all(abs(check_steady_state([Params.beta(end)*100; steadystate.Y])) < 1e-4))

[R, w, tau, dividend] = steadystateprices(steadystate.Y);


parstst = steadystatepolicy(steadystate.Y);
Pi = forwardmat(1,parstst,R,w,tau,dividend);
Dstst = invdistr(Pi);
clear Pi;


steadystate.w = w;
steadystate.R = R;
steadystate.tau = tau;
steadystate.dividend = dividend;

%% ZLB

T = 250;

%initial guess
Y = steadystate.Y;
wage = 1/Params.mu;
N = Y;
C = Y;
dividend = Y - wage*N;
S = 1;
ii = Params.Rbar-1;
ppi = 1;
pstar = 1;
pbarB = Y/(1 - Params.beta * (1-Params.theta));
pbarA = pbarB;
laborwedge = Params.psi1*N^(Params.psi2)/(C^(-Params.sigma)*wage);
eulerwedge = 1/(Params.beta *(1+ii));

X0 = repmat([Y; C; N; ppi; ii; S; dividend; pstar; pbarA; pbarB; wage],1,T);

    
%%
for iextended = 0:1
    
    X = X0; 
       
    %exog var
    extended = logical(iextended);
    disp(['------- extended = ' num2str(extended) ' --------'])
    
    betashock = ones(1,T);
    mpshock = zeros(1,T);
    
   
    
    betashock(2:34) = 1.001638;  % standard case
    %betashock(2:34) = 1.00175;  % high-risk case
    if extended
        mpshock(2:24) = -100;
        mpshock(25) = -0.00092;
       
        % no output gap at t = 0  
        %mpshock(25) = -100;
        %mpshock(26) = -0.0014;
       
        
%         mpshock(2:12) = -100;
%         mpshock(14:17) = -100;
    end
    
%     %SmallerLongerShock
%     betashock(2:33) = 1.0015;
%     if extended
%         mpshock(2:22) = -100;
%     end


    if numel(laborwedge) == 1
        laborwedge = repmat(laborwedge,1,T);
        eulerwedge = repmat(eulerwedge,1,T);
    end
    eps = [mpshock; betashock;laborwedge;eulerwedge];
    
    DetSimul_maxit = 100;
    
    Dampening = 0.25;
    DampeningThresh = 0.001;
    
    for it = 1:30
        
        
        %solve for prices
        disp('Updating prices')
        X  = DetSimul( X, eps, @ZLB_equ, DetSimul_maxit, Dampening, DampeningThresh);
        
        % update household behavior -> wedges
        ii = X(5,:);
        wage = X(11,:);
        C = X(2,:);
        N = X(3,:);
        ppi = X(4,:);
        dividend = X(7,:);
        [update, laborwedge, eulerwedge] = transition_hhld_update(laborwedge, eulerwedge, ii, wage, C, N, ppi, dividend, betashock, parstst, Dstst);
        
        eps(3:4,:) = [laborwedge; eulerwedge];
        
        if update < 5e-6,
            break
        end
    end
    
    Y = X(1,:);
    ii = X(5,:);
    wage = X(11,:);
    ppi = X(4,:);
    
    if extended
        fname = 'Results/transition_ZLB_extended.mat';
    else
        fname = 'Results/transition_ZLB_naive.mat';
    end
    save(fname,'Y','ii','wage','ppi');
    
end