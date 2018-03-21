function [residual, g1, g2, g3] = Q2_dynamic(y, x, params, steady_state, it_)
%
% Status : Computes dynamic model for Dynare
%
% Inputs :
%   y         [#dynamic variables by 1] double    vector of endogenous variables in the order stored
%                                                 in M_.lead_lag_incidence; see the Manual
%   x         [nperiods by M_.exo_nbr] double     matrix of exogenous variables (in declaration order)
%                                                 for all simulation periods
%   steady_state  [M_.endo_nbr by 1] double       vector of steady state values
%   params    [M_.param_nbr by 1] double          vector of parameter values in declaration order
%   it_       scalar double                       time period for exogenous variables for which to evaluate the model
%
% Outputs:
%   residual  [M_.endo_nbr by 1] double    vector of residuals of the dynamic model equations in order of 
%                                          declaration of the equations.
%                                          Dynare may prepend auxiliary equations, see M_.aux_vars
%   g1        [M_.endo_nbr by #dynamic variables] double    Jacobian matrix of the dynamic model equations;
%                                                           rows: equations in order of declaration
%                                                           columns: variables in order stored in M_.lead_lag_incidence followed by the ones in M_.exo_names
%   g2        [M_.endo_nbr by (#dynamic variables)^2] double   Hessian matrix of the dynamic model equations;
%                                                              rows: equations in order of declaration
%                                                              columns: variables in order stored in M_.lead_lag_incidence followed by the ones in M_.exo_names
%   g3        [M_.endo_nbr by (#dynamic variables)^3] double   Third order derivative matrix of the dynamic model equations;
%                                                              rows: equations in order of declaration
%                                                              columns: variables in order stored in M_.lead_lag_incidence followed by the ones in M_.exo_names
%
%
% Warning : this file is generated automatically by Dynare
%           from model file (.mod)

%
% Model equations
%

residual = zeros(7, 1);
T10 = exp(y(4))^(-1);
T15 = params(2)*exp(y(10))^(-1);
T24 = params(1)*exp(y(12))*exp(y(9))^(params(1)-1);
T28 = exp(y(11))^(1-params(1));
T29 = T24*T28;
T32 = 1+T29-params(4);
T41 = exp(y(2))^params(1);
T46 = exp(y(6))^(-params(1));
T63 = exp(y(6))^(1-params(1));
T64 = exp(y(8))*T41*T63;
lhs =T10;
rhs =T15*T32;
residual(1)= lhs-rhs;
lhs =0;
rhs =T10*(1-params(1))*exp(y(8))*T41*T46-params(3)*exp(y(6))^(1/params(7));
residual(2)= lhs-rhs;
lhs =exp(y(9));
rhs =exp(y(2))*(1-params(4))+exp(y(5));
residual(3)= lhs-rhs;
lhs =exp(y(4))+exp(y(5));
rhs =T64;
residual(4)= lhs-rhs;
lhs =y(8);
rhs =params(5)*y(1)+x(it_, 1);
residual(5)= lhs-rhs;
lhs =exp(y(3));
rhs =T64;
residual(6)= lhs-rhs;
lhs =exp(y(7));
rhs =T46*T41*(1-params(1))*exp(y(8));
residual(7)= lhs-rhs;
if nargout >= 2,
  g1 = zeros(7, 13);

  %
  % Jacobian matrix
  %

T82 = exp(y(4))*getPowerDeriv(exp(y(4)),(-1),1);
T95 = exp(y(6))*getPowerDeriv(exp(y(6)),(-params(1)),1);
T105 = (-(exp(y(8))*T41*exp(y(6))*getPowerDeriv(exp(y(6)),1-params(1),1)));
T120 = exp(y(2))*getPowerDeriv(exp(y(2)),params(1),1);
  g1(1,4)=T82;
  g1(1,10)=(-(T32*params(2)*exp(y(10))*getPowerDeriv(exp(y(10)),(-1),1)));
  g1(1,11)=(-(T15*T24*exp(y(11))*getPowerDeriv(exp(y(11)),1-params(1),1)));
  g1(1,12)=(-(T15*T29));
  g1(1,9)=(-(T15*T28*params(1)*exp(y(12))*exp(y(9))*getPowerDeriv(exp(y(9)),params(1)-1,1)));
  g1(2,4)=(-(T46*T41*exp(y(8))*(1-params(1))*T82));
  g1(2,6)=(-(T10*(1-params(1))*exp(y(8))*T41*T95-params(3)*exp(y(6))*getPowerDeriv(exp(y(6)),1/params(7),1)));
  g1(2,8)=(-(T10*(1-params(1))*exp(y(8))*T41*T46));
  g1(2,2)=(-(T46*T10*(1-params(1))*exp(y(8))*T120));
  g1(3,5)=(-exp(y(5)));
  g1(3,2)=(-(exp(y(2))*(1-params(4))));
  g1(3,9)=exp(y(9));
  g1(4,4)=exp(y(4));
  g1(4,5)=exp(y(5));
  g1(4,6)=T105;
  g1(4,8)=(-T64);
  g1(4,2)=(-(T63*exp(y(8))*T120));
  g1(5,1)=(-params(5));
  g1(5,8)=1;
  g1(5,13)=(-1);
  g1(6,3)=exp(y(3));
  g1(6,6)=T105;
  g1(6,8)=(-T64);
  g1(6,2)=(-(T63*exp(y(8))*T120));
  g1(7,6)=(-(T41*(1-params(1))*exp(y(8))*T95));
  g1(7,7)=exp(y(7));
  g1(7,8)=(-(T46*T41*(1-params(1))*exp(y(8))));
  g1(7,2)=(-(T46*(1-params(1))*exp(y(8))*T120));

if nargout >= 3,
  %
  % Hessian matrix
  %

  g2 = sparse([],[],[],7,169);
if nargout >= 4,
  %
  % Third order derivatives
  %

  g3 = sparse([],[],[],7,2197);
end
end
end
end
