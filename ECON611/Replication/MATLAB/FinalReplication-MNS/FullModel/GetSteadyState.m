function [steadystate, parstst, Dstst] = GetSteadyState()

global Params;

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


end