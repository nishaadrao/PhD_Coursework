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

names = {'R','A', 'S', 'wage','dividend','ppi','Y','C','N','L','pbarB','pbarA','pstar'};

nvar = length(names);
stst = zeros(nvar,1);

for i_ = 1:nvar
    eval(['stst(i_) = ' names{i_} ';'])
    eval(['ind_' names{i_} ' = i_;'])
end

clear(names{:});

%% Solve for transition


T = 200;

%Initial guesses
[R, w] = steadystateprices();
wagepath = repmat(w,1,T);
dividendpath = repmat(stst(ind_dividend),1,T);

Spath = ones(1,T);


%policy
Rpath =  repmat(R,1,T);
HORIZ = 20;
Rpath(HORIZ+2) = 1;
%Rpath(81) = 1.01;



eqm = SolveForTransitionCompMkts( Rpath ,  wagepath, dividendpath, Spath, stst, names );



figOutput = figure;
CMOutputGE = eqm(ind_Y,2:T-1)./stst(ind_Y)-1;
plot(0:T-3,CMOutputGE)
title('equilibrium output under complete markets')
saveas(figOutput,'Figures/CompleteMarketsOutput.pdf')
save('Results/CMOutputGE.mat','CMOutputGE');


figure;
CMInflationGE = eqm(ind_ppi,2:T-1)-1;
plot(0:T-3,CMInflationGE)
title('equilibrium output under complete markets')
saveas(gcf,'Figures/CompleteMarketsOutput.pdf')
save('Results/CMInflationGE.mat','CMInflationGE');

% key moments
disp('change in output and inflation at date 0')
disp([eqm(ind_Y,2)/stst(ind_Y)  eqm(ind_ppi,2)]*10000-10000)
disp('change in output and inflation at date 20')
disp([eqm(ind_Y,22)/stst(ind_Y)  eqm(ind_ppi,22)]*10000-10000)

%% Horizon figure

T = 200;


%load pricepath;
wagepath = repmat(stst(ind_wage),1,T);
dividendpath = repmat(stst(ind_dividend),1,T);
Spath = ones(1,T);


% list of horizons at which change occurs
horiz = [1:30]';


HorizonTable = zeros(length(horiz),3);

for ih = 1:length(horiz)

    Rpath =  repmat(stst(ind_R),1,T);
    Rpath(horiz(ih)+2) = 1;
    %Rpath(horiz(ih)+1) = 1.01;
    eqm  = SolveForTransitionCompMkts( Rpath ,  wagepath, dividendpath, Spath, stst, names);
    HorizonTable(ih,:) = [horiz(ih) eqm(ind_Y,2)/stst(ind_Y) eqm(ind_ppi,2)];
    
    
end
save('Results/HorizonTable_CM_Working.mat','HorizonTable');





