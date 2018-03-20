function [residual, g1, g2, g3] = Q2_static(y, x, params)
%
% Status : Computes static model for Dynare
%
% Inputs : 
%   y         [M_.endo_nbr by 1] double    vector of endogenous variables in declaration order
%   x         [M_.exo_nbr by 1] double     vector of exogenous variables in declaration order
%   params    [M_.param_nbr by 1] double   vector of parameter values in declaration order
%
% Outputs:
%   residual  [M_.endo_nbr by 1] double    vector of residuals of the static model equations 
%                                          in order of declaration of the equations.
%                                          Dynare may prepend or append auxiliary equations, see M_.aux_vars
%   g1        [M_.endo_nbr by M_.endo_nbr] double    Jacobian matrix of the static model equations;
%                                                       columns: variables in declaration order
%                                                       rows: equations in order of declaration
%   g2        [M_.endo_nbr by (M_.endo_nbr)^2] double   Hessian matrix of the static model equations;
%                                                       columns: variables in declaration order
%                                                       rows: equations in order of declaration
%   g3        [M_.endo_nbr by (M_.endo_nbr)^3] double   Third derivatives matrix of the static model equations;
%                                                       columns: variables in declaration order
%                                                       rows: equations in order of declaration
%
%
% Warning : this file is generated automatically by Dynare
%           from model file (.mod)

residual = zeros( 6, 1);

%
% Model equations
%

T10 = exp(y(1))^(-1);
T21 = params(1)*exp(y(5))*exp(y(2))^(params(1)-1);
T25 = exp(y(4))^(1-params(1));
T29 = 1+T21*T25-params(4);
T34 = exp(y(2))^params(1);
T37 = exp(y(4))^(-params(1));
lhs =T10;
rhs =T10*params(2)*T29;
residual(1)= lhs-rhs;
lhs =0;
rhs =exp(y(5))*T10*(1-params(1))*T34*T37-exp(y(4))*params(3);
residual(2)= lhs-rhs;
lhs =exp(y(2));
rhs =exp(y(2))*(1-params(4))+exp(y(3));
residual(3)= lhs-rhs;
lhs =exp(y(1))+exp(y(3));
rhs =T25*exp(y(5))*T34;
residual(4)= lhs-rhs;
lhs =y(5);
rhs =y(5)*params(5)+x(1);
residual(5)= lhs-rhs;
lhs =exp(y(6));
rhs =T25*exp(y(5))*T34;
residual(6)= lhs-rhs;
if ~isreal(residual)
  residual = real(residual)+imag(residual).^2;
end
if nargout >= 2,
  g1 = zeros(6, 6);

  %
  % Jacobian matrix
  %

T62 = exp(y(1))*getPowerDeriv(exp(y(1)),(-1),1);
T78 = exp(y(2))*getPowerDeriv(exp(y(2)),params(1),1);
T88 = exp(y(4))*getPowerDeriv(exp(y(4)),1-params(1),1);
  g1(1,1)=T62-T29*params(2)*T62;
  g1(1,2)=(-(T10*params(2)*T25*params(1)*exp(y(5))*exp(y(2))*getPowerDeriv(exp(y(2)),params(1)-1,1)));
  g1(1,4)=(-(T10*params(2)*T21*T88));
  g1(1,5)=(-(T10*params(2)*T21*T25));
  g1(2,1)=(-(T37*T34*exp(y(5))*(1-params(1))*T62));
  g1(2,2)=(-(T37*exp(y(5))*T10*(1-params(1))*T78));
  g1(2,4)=(-(exp(y(5))*T10*(1-params(1))*T34*exp(y(4))*getPowerDeriv(exp(y(4)),(-params(1)),1)-exp(y(4))*params(3)));
  g1(2,5)=(-(exp(y(5))*T10*(1-params(1))*T34*T37));
  g1(3,2)=exp(y(2))-exp(y(2))*(1-params(4));
  g1(3,3)=(-exp(y(3)));
  g1(4,1)=exp(y(1));
  g1(4,2)=(-(T25*exp(y(5))*T78));
  g1(4,3)=exp(y(3));
  g1(4,4)=(-(exp(y(5))*T34*T88));
  g1(4,5)=(-(T25*exp(y(5))*T34));
  g1(5,5)=1-params(5);
  g1(6,2)=(-(T25*exp(y(5))*T78));
  g1(6,4)=(-(exp(y(5))*T34*T88));
  g1(6,5)=(-(T25*exp(y(5))*T34));
  g1(6,6)=exp(y(6));
  if ~isreal(g1)
    g1 = real(g1)+2*imag(g1);
  end
if nargout >= 3,
  %
  % Hessian matrix
  %

  g2 = sparse([],[],[],6,36);
if nargout >= 4,
  %
  % Third order derivatives
  %

  g3 = sparse([],[],[],6,216);
end
end
end
end
