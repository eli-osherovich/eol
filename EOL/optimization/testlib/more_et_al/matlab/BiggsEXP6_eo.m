function [fvec, Jac, hessV, hessVfunc, x0_std] = BiggsEXP6_eo(x, ~, ...
                                                  M)
%BIGGSEXP6 - Biggs EXP6 function.

% Number of variables: N = 6
% Number of equations: M = any M>=N (M=13 be default)
% Standard starting point x0_std = [1; 2; 1; 1; 1; 1]
% Min function value = 5.65565E-3 (if M=13)
% Argmin x =
% Additional extrema:
% f=0 at x=[1; 10; 1; 5; 4; 3]


% Copyright 2010 Eli Osherovich.


% standard starting point
x0_std = [1; 2; 1; 1; 1; 1];

if 5 == nargout
    % return empty values if x0_std requested
    fvec = [];
    Jac = [];
    hessV = [];
    hessVfunc = [];
    return;
end

% set M if not provided
if 3 > nargin
    M = 13;
end

% initialize some constants
t = (1:M)'/10;
y = exp(-t) - 5*exp(-10*t) + 3*exp(-4*t);


fvec = x(3)*exp(-t*x(1)) - x(4)*exp(-t*x(2)) + x(6)*exp(-t*x(5)) - ...
       y;



if nargout > 1 % Jacobian requested
    Jac = zeros(M, 6);
    Jac(:, 1) = -x(3)*t.*exp(-t*x(1));
    Jac(:, 2) = x(4)*t.*exp(-t*x(2));
    Jac(:, 3) = exp(-t*x(1));
    Jac(:, 4) = -exp(-t*x(2));
    Jac(:, 5) = -x(6)*t.*exp(-t*x(5));
    Jac(:, 6) = exp(-t*x(5));
    
    switch nargout
      case {3, 4} % Hessian or Hessian-vector product
        error('EOL:BiggsEXP6:HessianNotSupported', ['Hessian or ' ...
                            'Hessian-vector product are not supported']);
      case 5
        hessV = [];
        hessVfunc = [];
    end
end
