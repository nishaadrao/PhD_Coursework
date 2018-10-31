function FullModel(DoHorizon)

if ~exist('DoHorizon','var')
    DoHorizon = false;
end

setpath();
run('FullModel/parameters/initparams');


%% Adjust the discount factor to target a certain level of assets


Y0 = 0.6;
optset('broyden','tol',1e-6);
x = broyden(@check_steady_state,[100*Params.beta(end);Y0]);

if Params.beta_heterog
    Params.beta = kron(x(1)/100+[-Params.betad 0],ones(1,Params.npp/Params.nbeta));
else
    Params.beta = x(1)/100;
end
steadystate.Y = x(2);
Params.B = steadystate.Y * Params.asset_target;

assert(all(abs(check_steady_state([Params.beta(end)*100; steadystate.Y])) < 1e-3))

[R, w, tau, dividend] = steadystateprices(steadystate.Y);


parstst = steadystatepolicy(steadystate.Y);
Pi = forwardmat(1,parstst,R,w,tau,dividend);
Dstst = invdistr(Pi);
clear Pi;


steadystate.w = w;
steadystate.R = R;
steadystate.tau = tau;
steadystate.dividend = dividend;



%% Analyze the steady state

if exist('Figures/SteadyStateProperties.txt','file')
    delete('Figures/SteadyStateProperties.txt')
end
diary('Figures/SteadyStateProperties.txt')

b = Params.knotDistrK + Params.borrowCon;
CDF = cumsum(sum(reshape(Dstst,Params.ndstst,Params.npp),2));
[~,i50] = min(abs(CDF-0.5));




disp(' ');
GDP_PER_CAPITA_2010 = 48377.39;
disp(['Median assets as frac of annual income = ' num2str(b(i50)/(steadystate.Y*4))])
disp(['Median assets in 2010 $ = ' num2str(b(i50)/(steadystate.Y*4)*GDP_PER_CAPITA_2010)])

disp(' ');
disp('Fraction of total population with low wealth');
disp('income type     share of pop.       fraction with low weatlh ');
lowWealth = max(Params.endow);
for ip = 1:Params.npp
    [lwshare, popshare] =expect_low_wealth(Dstst,ip,lowWealth);
    disp([num2str(ip) '           ' num2str(popshare) '             ' num2str(lwshare)])
end

disp(' ');
aggassets = expect_k(Dstst);
disp(['aggregate assets to annual GDP ' num2str( aggassets/(steadystate.Y*4))]);
disp(['aggregate debt to annual GDP ' num2str(expect_debt(Dstst)/(steadystate.Y*4))]);
disp(['Fraction with positive assets ' num2str(expect_debt(Dstst,true))]);


WlthFig.x = (Params.knotDistrK+Params.borrowCon)/steadystate.Y/4;
WlthFig.y = cumsum(sum(reshape(Dstst,1000,Params.npp),2));
save(['Figures/WlthDistribution_assets_' num2str(Params.asset_target) '.mat'],'WlthFig');



diary off

clearvars -except Dstst parstst Params steadystate DoHorizon;




%% Compute transition in GE
T = 200;

%Initial guesses
[R, w] = steadystateprices();


%load pricepath;
wagepath = repmat(w,1,T);
dividendpath = repmat(steadystate.dividend,1,T);
Spath = ones(1,T);


%policy
Rpath =  repmat(R,1,T);
HORIZ = 20;
Rpath(HORIZ+2) = 1;
eqm  = SolveForTransition( Rpath, wagepath, dividendpath, Spath, steadystate, parstst, Dstst);


figy = figure;
IMOutputGE = eqm.Y(2:end-1)/steadystate.Y-1;
plot([0:T-3],100*IMOutputGE)
xlim([0 40])
xlabel('Quarter')
ylabel('Percentage points')
saveas(figy,'Figures/output.pdf')
save('Results/IMOutputGE.mat','IMOutputGE')

figpi = figure;
IMInflationGE = eqm.ppi(2:end-1)-1;
plot([0:T-3],IMInflationGE)
title('equilibrium inflation')
saveas(figy,'Figures/inflation.pdf')
save('Results/IMInflationGE.mat','IMInflationGE')

figwealth = figure;
WealthShares = eqm.WealthShares(:,2:40)'./(ones(39,1)*eqm.WealthShares(:,1)');
plot(WealthShares);
legend('low prod.', 'med. prod.', 'high prod.')
title('evolution of wealth shares')
saveas(figwealth, 'Figures/wealthshares.pdf');
save('Results/WealthSharesGE.mat','WealthShares');

% key moments
disp('change in output and inflation at date 0')
disp([eqm.Y(2)/steadystate.Y  eqm.ppi(2)]*10000-10000)
disp('change in output and inflation at date 20')
disp([eqm.Y(22)/steadystate.Y  eqm.ppi(22)]*10000-10000)




%% Create data for plots over different horizons


if DoHorizon
    T = 200;
    
    
    %load pricepath;
    wagepath = repmat(steadystate.w,1,T);
    dividendpath = repmat(steadystate.dividend,1,T);
    Spath = ones(1,T);
    
    
    % list of horizons at which change occurs
    horiz = [1:2:41]';
    
    
    HorizonTable = zeros(length(horiz),3);
    
    for ih = 1:length(horiz)
        tic
        Rpath =  repmat(steadystate.R,1,T);
        Rpath(horiz(ih)+2) = 1;
        %Rpath(horiz(ih)+1) = 1.01;
        eqm  = SolveForTransition( Rpath, wagepath, dividendpath, Spath, steadystate, parstst, Dstst);
        HorizonTable(ih,:) = [horiz(ih) eqm.Y(2)/steadystate.Y eqm.ppi(2)];
        save('Results/HorizonTable_IM_Working.mat','HorizonTable');
        toc
        
    end
    
    
end


