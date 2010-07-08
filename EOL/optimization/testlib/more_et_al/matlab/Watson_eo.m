function [fvec, Jac, hessV, hessVfunc, x0_std] = Watson_eo(x, ~, ~, N)
%Watson - Watson function.

% Number of variables: N = any (2<= N <= 31) (N=12 by default)
% Number of equations: M = 31
% Standard starting point x0_std = [0;0;...;0]
% Min function value = 
% Argmin x =
% Additional extrema:
% f=2.28767E-3 if N=6
% f=1.39976E-6 if N=9
% f=4.72238E-10 if N=12 


% Copyright 2010 Eli Osherovich.

% set N
if 4 > nargin 
    N = 12;
end

% set M
M = 31;

% standard starting point
x0_std = zeros(N,1);

if 5 == nargout
    % return empty values if x0_std requested
    fvec = [];
    Jac = [];
    hessV = [];
    hessVfunc = [];
    return;
end

% initialize constants
t = (1:29)'/29;
J = (1:N)';


% compute f values
fvec = zeros(M,1);
for i = 1:29
    fvec(i) = sum((J-1).*x.*t(i).^(J-2)) - sum(x.*t(i).^(J-1))^2 - 1;
end
fvec(30) = x(1);
fvec(31) = x(2) - x(1)^2 - 1;

if nargout > 1 % Jacobian requested
    Jac = zeros(M,N);
    for i=1:29
        Jac(i,:) = (J-1).*t(i).^(J-2) - 2*sum(x.*t(i).^(J-1))*t(i).^(J-1);
    end
    Jac(30, 1) = 1;
    Jac(31, 1) = -2*x(1);
    Jac(31, 2) = 1;
    
    switch nargout
      case {3, 4} % Hessian or Hessian-vector product
        error('EOL:Watson:HessianNotSupported', ['Hessian or Hessian-vector ' ...
                            'product are not supported']);
      case 5
        hessV = [];
        hessVfunc = [];
    end
end
