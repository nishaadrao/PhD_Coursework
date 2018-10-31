close all
clear all
clc

run('FullModel/parameters/initparams');

% adjust psi1 for rep agent.  Then rep agent supplies labor with skill
% normalized to one.
Params.psi1 = Params.psi1 * dot(invdistr(Params.Trans'), Params.endow.^(1+1/Params.psi2));
Params.beta = 1/Params.Rbar;


%% Solve for steady state



A = 1;
ppi = 1;
wage = 1/Params.mu;
Y = (wage/Params.psi1)^(1/(Params.sigma + Params.psi2));
R = 1/Params.beta;
C = Y;
S = 1;
N = S*Y/A;
L = N;
pbarB = Y/(1- Params.beta * (1-Params.theta));
pbarA = wage*Params.mu*Y/(1-Params.beta * (1-Params.theta));
pstar = 1;
dividend = Y - wage*N;
ii = Params.Rbar-1;


%% ZLB

T = 250;

X0 = repmat([Y; C; N; ppi; ii; S; dividend; pstar; pbarA; pbarB; wage],1,T);

%%
for iextended = 0:1
    
    X = X0;
    
    %exog var
    extended = logical(iextended);
    disp(['------- extended = ' num2str(extended) ' --------'])
    
    betashock = ones(1,T);
    mpshock = zeros(1,T);
    
    %baseline
    betashock(2:34) = 1.0014825;
    %betashock(2:21) = 1.00238;  %high assets
    %betashock(2:21) = 1.0018;
    %betashock(2:21) = 1.0022;
    if extended
        mpshock(2:24) = -100;
        mpshock(25) = -0.00092;
        
        %mpshock(2:12) = -100;
        %mpshock(14:17) = -100;
    end

%     %SmallerLongerShock
%     betashock(2:33) = 1.001367;
%     if extended
%         mpshock(2:22) = -100;
%     end

    
    eps = [mpshock; betashock];
    
    DetSimul_maxit = 60;
    %solve for transition
    X  = DetSimul( X, eps, @ZLB_cmpMkts_equ, DetSimul_maxit,0.25,0.001);
    
    Y = X(1,:);
    ii = X(5,:);
    wage = X(11,:);
    ppi = X(4,:);
    
    if extended
        fname = 'Results/transition_ZLB_cmp_mkts_extended.mat';
    else
        fname = 'Results/transition_ZLB_cmp_mkts_naive.mat';
    end
    save(fname,'Y','ii','wage','ppi');
    
end