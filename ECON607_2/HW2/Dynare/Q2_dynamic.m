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

residual = zeros(6, 1);
T10 = exp(y(3))^(-1);
T15 = params(2)*exp(y(9))^(-1);
T24 = params(1)*exp(y(11))*exp(y(4))^(params(1)-1);
T28 = exp(y(10))^(1-params(1));
T32 = 1+T24*T28-params(4);
T41 = exp(y(1))^params(1);
T46 = exp(y(6))^(-params(1));
T60 = exp(y(6))^(1-params(1));
T61 = exp(y(7))*T41*T60;
lhs =T10;
rhs =T15*T32;
residual(1)= lhs-rhs;
lhs =0;
rhs =T10*(1-params(1))*exp(y(7))*T41*T46-exp(y(6))*params(3);
residual(2)= lhs-rhs;
lhs =exp(y(4));
rhs =exp(y(1))*(1-params(4))+exp(y(5));
residual(3)= lhs-rhs;
lhs =exp(y(3))+exp(y(5));
rhs =T61;
residual(4)= lhs-rhs;
lhs =y(7);
rhs =params(5)*y(2)+x(it_, 1);
residual(5)= lhs-rhs;
lhs =exp(y(8));
rhs =T61;
residual(6)= lhs-rhs;
if nargout >= 2,
  g1 = zeros(6, 12);

  %
  % Jacobian matrix
  %

T73 = exp(y(3))*getPowerDeriv(exp(y(3)),(-1),1);
T85 = exp(y(1))*getPowerDeriv(exp(y(1)),params(1),1);
T108 = (-(exp(y(7))*T41*exp(y(6))*getPowerDeriv(exp(y(6)),1-params(1),1)));
  g1(1,3)=T73;
  g1(1,9)=(-(T32*params(2)*exp(y(9))*getPowerDeriv(exp(y(9)),(-1),1)));
  g1(1,4)=(-(T15*T28*params(1)*exp(y(11))*exp(y(4))*getPowerDeriv(exp(y(4)),params(1)-1,1)));
  g1(1,10)=(-(T15*T24*exp(y(10))*getPowerDeriv(exp(y(10)),1-params(1),1)));
  g1(1,11)=(-(T15*T24*T28));
  g1(2,3)=(-(T46*T41*exp(y(7))*(1-params(1))*T73));
  g1(2,1)=(-(T46*T10*(1-params(1))*exp(y(7))*T85));
  g1(2,6)=(-(T10*(1-params(1))*exp(y(7))*T41*exp(y(6))*getPowerDeriv(exp(y(6)),(-params(1)),1)-exp(y(6))*params(3)));
  g1(2,7)=(-(T10*(1-params(1))*exp(y(7))*T41*T46));
  g1(3,1)=(-(exp(y(1))*(1-params(4))));
  g1(3,4)=exp(y(4));
  g1(3,5)=(-exp(y(5)));
  g1(4,3)=exp(y(3));
  g1(4,1)=(-(T60*exp(y(7))*T85));
  g1(4,5)=exp(y(5));
  g1(4,6)=T108;
  g1(4,7)=(-T61);
  g1(5,2)=(-params(5));
  g1(5,7)=1;
  g1(5,12)=(-1);
  g1(6,1)=(-(T60*exp(y(7))*T85));
  g1(6,6)=T108;
  g1(6,7)=(-T61);
  g1(6,8)=exp(y(8));

if nargout >= 3,
  %
  % Hessian matrix
  %

  g2 = sparse([],[],[],6,144);
if nargout >= 4,
  %
  % Third order derivatives
  %

  g3 = sparse([],[],[],6,1728);
end
end
end
end
