// This is a Dynare file for a basic New Keynesian Model as discussed in 
// Jordi Galí's textbook "Monetary Policy, Inflation, and the Business Cycle" (2008).
//
// The model is simulated in order to see the evolution of the deviations from steady state for endogenous variables.
// The driving force is a technology shock. Monetary policy shocks can be considered once the shock e_v is turned on.
//
// All equations are expressed in their logarithmic form. Results are to be intepreted as deviations from steady state.
//
// Enjoy.
//
// Author:
// Johannes Fritz, University of St.Gallen (Switzerland), johannes.fritz[A]unisg.ch.

var y, pi, i, a, rn, n, m, v;
varexo e_a e_v;

parameters alpha, beta, theta, sigma, phi, rho, phi_pi, phi_y, rho_a, rho_v, lambda, kappa, psi, epsilon, eta;

// These values follow the (corrected) estimation of Galí, Gertler & López-Salido (2001/3) for the US, specification (2) in Table 1.
alpha=.27;              // If =0 we have a CRS production technology. Else it's decreasing returns to scale (see model equation 5).
epsilon=1.5;            // Elasticity of substitution derived from the markup forumla m=log(epsilon/(epsilon-1)). Using m=1.1.
beta=0.923;             // The discount factor. 
theta=0.698;            // Measure of price stickiness. If =0 then prices are flexible.
lambda=0.154;           // or alternatively derived endogneously through lambda=(theta^(-1))*(1-theta)*(1-beta*theta)*(1-alpha)/(1-alpha+alpha*epsilon).
rho=-log(beta);         // Real interest rate in the steady state (no shocks).

sigma=1;                // Coefficient of risk aversion.
phi=1;                  // Elasticity of labor supply.
phi_pi=1.5;             // Sensitivity of the central bank with respect to inflation.
phi_y=0.5/4;              // Sensitivity of the central bank with respect to the output gap.
rho_a=0.975;            // Persistence of the technology shock.
rho_v=0.5;              // Persistence of the monetary policy shock.
eta=4;                  // Elasticity of the money demand with respect to the nominal interest rate (see Eq.6).


// The next two parameters are generated for the solution of the model. Note that when alpha=0, these equations get much easier.
kappa=lambda*(sigma+(phi+alpha)/(1-alpha));
psi=(1+phi)*((sigma+phi+alpha*(1-sigma))^(-1));


model;                              // All equations are stated in log form.
y=y(+1)-1/sigma*(i-pi(+1)-rn);      // Eq. 1: The Dynamic IS equation.
pi=beta*pi(+1)+kappa*y;             // Eq. 2: The New Keynesian Philips Curve.
rn=rho+sigma*psi*(rho_a-1)*a;       // Eq. 3: The evolution of the natural rate of interest.
i=rho+phi_pi*pi+phi_y*y+v;          // Eq. 4: The interest rate rule of the central bank.
y=a+(1-alpha)*n;                    // Eq. 5: The production function consisting of technology and labor. This relationship is only true up to a 1st order approximation.
m=pi+y-eta*(i); 	            // Eq. 6: Ad-hoc money demand.
a=rho_a*a(-1)+e_a;                  // Eq. 7: Technology shocks follow an AR(1) process with persistence rho_a.
v=rho_v*v(-1)+e_v;                  // Eq. 8: Monetary policy follow an AR(1) process with persistence rho_v.
end;

initval;
y=0;
m=0;
n=0;
pi=0;
i=rho;
rn=rho;
a=0;
v=0;
e_a=0;
e_v=0;
end;

steady;
check;

shocks;
var e_a;
stderr 0.0072;
var e_v;
stderr 0;        // The monetary policy shock is turned off.
end;


// The above equations only hold up to a first order approximation. Thus order=1 for the simulation.
stoch_simul(irf=12, nofunctions, order=1) n y i rn pi a m v; 