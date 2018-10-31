
global Params;


%----- Options for alternative calibrations ----:
HIGH_RISK = false;
HIGH_ASSET = false;
%---------------------------------------------


Params.Rbar = 1.005;
Params.sigma = 2;     %risk aversion


Params.mu = 1.2;

Params.psi1 = 1;
Params.psi2 = 2;
%Params.psi2 = 1/0.9;

Params.borrowCon = 0; 
%Params.borrowCon = -5/12*4 *  0.4727;  % 5 times monthly labor income of low-skill households


Params.theta = 0.15;


Params.phiP = 1.5;

if HIGH_ASSET
    Params.asset_target = 15.14;  % high-asset
else
    Params.asset_target = 5.5;  % baseline
end
%Params.asset_target = 2.5;  %matches median liquid wealth






%% Household heterogeneity and income process




%p1 = 0.0759;
%e2 = 0.96;
%Params.wage_rho = 0.961797292688304;
%Params.wage_sigma = 0.177099384984742^2/(1-Params.wage_rho^2);
%[ Params.Trans, Params.endow ] = skillCalib( p1, e2, Params.wage_sigma, Params.wage_rho );

%rho_z = 0.98;

Params.npp = 3;
%[Params.Trans, Params.endow] = rouwen(rho_z, 0, sqrt(0.4), Params.npp);

if ~HIGH_RISK   %baseline
    uncondSTD = sqrt(0.01695/(1-0.96566^2));
    [Params.Trans, Params.endow] = rouwen(0.96566, 0, uncondSTD, Params.npp);
    %[Params.Trans, Params.endow] = rouwen(0.98, 0, uncondSTD, Params.npp);  %High persistence
elseif Params.asset_target < 10  % High-risk case
    [Params.Trans, Params.endow] = rouwen(0.96566, 0, sqrt(0.033/(1-0.96566^2)), Params.npp);
else % high-risk and high-asset case
    [Params.Trans, Params.endow] = rouwen(0.96566, 0, sqrt(0.024/(1-0.96566^2)), Params.npp);
end

Params.endow = exp(Params.endow)';
Params.Trans = Params.Trans';
Params.Trans_orig = Params.Trans;
Params.tax_weights = Params.endow;
%RICH PAY TAX:
Params.tax_weights(1:2) = 0;
%Params.tax_weights = [0.01 0.111 0.232] .* Params.endow;

Params.AvgTaxWeight = dot(invdistr(Params.Trans'), Params.tax_weights);



%% Preference heterog
Params.beta_heterog = false;

if Params.beta_heterog
    Params.nbeta = 3;
    betaHi = 0.994;
    Params.betad = [0.05 0.01];
    Params.beta = kron([betaHi-Params.betad betaHi], ones(1,Params.npp));
    beta_trans = [0.995  0.005 0; 0.0025 0.995 0.0025; 0 0.005 0.995];
    assert(all(abs(sum(beta_trans,2)-1.0)<1e-6))
    Params.Trans = kron(beta_trans,Params.Trans);
    Params.endow = repmat(Params.endow,1,Params.nbeta);
    Params.tax_weights = repmat(Params.tax_weights,1,Params.nbeta);
    Params.npp = Params.npp*Params.nbeta;
else
    Params.beta = 0.99;
end


%% Approximation of policy rules and asset positions

Params.nc = 200; %number of points in approximation of individual cons func


Params.ndstst = 1000; %number of points in wealth histogram


Params.bhmin = 0;  %min bond position for household
Params.bhmax = 75;  %max bond position for household





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% quadrature grid for_ wealth distribution:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Params.nStatesDistr = Params.ndstst*Params.npp; 
% number knot points:
ndk = Params.ndstst ;
[Params.knotDistrK,Params.logshift] = makeknotd(Params.bhmin,Params.bhmax,ndk);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% KNOT POINTS FOR_ CONSUMPTION POLYNOMIAL:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Params.xmin = 0.001;
Params.xmax = Params.bhmax;


Params.bgrid = logspaceshift(Params.xmin,Params.xmax,Params.nc,Params.logshift)';




